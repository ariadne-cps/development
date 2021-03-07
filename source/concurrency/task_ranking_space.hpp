/***************************************************************************
 *            concurrency/task_ranking_space.hpp
 *
 *  Copyright  2007-20  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file concurrency/task_raking_space.hpp
 *  \brief Class for handling ranking of results from parameter search of a task.
 */

#ifndef ARIADNE_TASK_RANKING_SPACE_HPP
#define ARIADNE_TASK_RANKING_SPACE_HPP

#include "../algebra/vector.hpp"
#include "../concurrency/task_ranking_parameter.hpp"
#include "../concurrency/task_execution_ranking.hpp"

namespace Ariadne {

using std::min, std::max;

template<class R> class TaskRankingSpace;
typedef ScoreType WeightType;
template<class R> using ParameterWeightsList = List<Pair<TaskRankingParameter<R>,WeightType>>;

template<class R>
class TaskRankingSpaceBuilder {
public:
    TaskRankingSpaceBuilder() = default;

    TaskRankingSpaceBuilder& add(TaskRankingParameter<R> const& parameter, WeightType const& weight = WeightType(1.0)) {
        ARIADNE_PRECONDITION(weight >= 0);
        _parameter_weights.push_back(Pair<TaskRankingParameter<R>,WeightType>(parameter, weight));
        return *this;
    }

    TaskRankingSpace<R> build() const;

  private:
    ParameterWeightsList<R> _parameter_weights;
};

template<class R>
class TaskRankingSpace : public WritableInterface {
    friend class TaskRankingSpaceBuilder<R>;
  public:
    typedef ScoreType WeightType;
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;

  protected:
    TaskRankingSpace(List<Pair<TaskRankingParameter<R>,WeightType>> const& parameter_weights)
        : _parameter_weights(parameter_weights) { }
  public:
    ParameterWeightsList<R> const& parameter_weights() const { return _parameter_weights; }

    TaskRankingSpace* clone() const { return new TaskRankingSpace(*this); }

    virtual OutputStream& _write(OutputStream& os) const { os << _parameter_weights; return os; }

    Bool has_critical_constraints() const {
        for (auto pw : _parameter_weights) if (pw.first.severity() == RankingConstraintSeverity::CRITICAL) return true;
        return false;
    }

    List<TaskRankingParameter<R>> failed_critical_constraints(InputType const& input, OutputType const& output) const {
        List<TaskRankingParameter<R>> result;
        for (auto pw : _parameter_weights) {
            if (pw.first.severity() == RankingConstraintSeverity::CRITICAL) {
                auto const& p = pw.first;
                auto rank = p.rank(input,output,DurationType(0));
                auto threshold = p.threshold(input,output,DurationType(0));
                if ((p.optimisation() == OptimisationCriterion::MINIMISE and rank > threshold) or
                    (p.optimisation() == OptimisationCriterion::MAXIMISE and rank < threshold)) {
                    result.push_back(p);
                }
            }
        }
        return result;
    }

    List<TaskExecutionRanking> rank(Map<ConfigurationSearchPoint,Pair<OutputType,DurationType>> const& data, InputType const& input) const {
        List<TaskExecutionRanking> result;

        ParameterWeightsList<R> kept_parameter_weights;
        for (auto pw : _parameter_weights) if (not pw.first.discard(input)) kept_parameter_weights.push_back(pw);

        // Compute the dimensions and initialise the min/max entries
        List<SizeType> dimensions;
        List<Pair<ScoreType,ScoreType>> scalar_min_max;
        List<Vector<Pair<ScoreType,ScoreType>>> vector_min_max;
        auto data_iter = data.cbegin();
        for (auto pw : kept_parameter_weights) {
            auto p = pw.first;
            auto dim = p.dimension(input);
            dimensions.push_back(dim);
            if (dim == 1) {
                auto val = p.rank(input, data_iter->second.first, data_iter->second.second);
                scalar_min_max.push_back(Pair<ScoreType,ScoreType>{val, val});
            } else {
                Vector<Pair<ScoreType,ScoreType>> vals(dim);
                for (SizeType i=0; i<dim; ++i) {
                    auto val = p.rank(input, data_iter->second.first, data_iter->second.second, i);
                    vals[i] = {val,val};
                }
                vector_min_max.push_back(vals);
            }
        }

        // Update the min/max on the remaining data entries
        ++data_iter;
        while (data_iter != data.cend()) {
            for (SizeType pw_idx=0; pw_idx < kept_parameter_weights.size(); ++pw_idx) {
                auto p = kept_parameter_weights[pw_idx].first;
                auto dim = dimensions[pw_idx];
                if (dim == 1) {
                    auto val = p.rank(input, data_iter->second.first, data_iter->second.second);
                    scalar_min_max[pw_idx] = {min(scalar_min_max[pw_idx].first,val), max(scalar_min_max[pw_idx].second,val)};
                } else {
                    for (SizeType i=0; i<dim; ++i) {
                        auto val = p.rank(input, data_iter->second.first, data_iter->second.second, i);
                        vector_min_max[pw_idx][i] = {min(vector_min_max[pw_idx][i].first,val),max(vector_min_max[pw_idx][i].second,val)};
                    }
                }
            }
            ++data_iter;
        }

        // Compute the score
        for (auto entry : data) {
            // If max != min, compute (val-min)/(max-min), in the vector case also divide by the dimension
            // sum/subtract each weighted cost based on the minimise/maximise objective
            ScoreType score(0);
            SizeType low_errors(0), high_errors(0);
            for (SizeType pw_idx=0; pw_idx < kept_parameter_weights.size(); ++pw_idx) {
                auto p = kept_parameter_weights[pw_idx].first;
                auto weight = kept_parameter_weights[pw_idx].second;
                auto dim = dimensions[pw_idx];
                ScoreType local_score(0);
                if (dim == 1) {
                    auto max_min_diff = scalar_min_max[pw_idx].second - scalar_min_max[pw_idx].first;
                    auto rank = p.rank(input, entry.second.first, entry.second.second);
                    if (p.uses_objective()) {
                        auto threshold = p.threshold(input, entry.second.first, entry.second.second);
                        if ((p.optimisation() == OptimisationCriterion::MINIMISE and rank > threshold) or
                            (p.optimisation() == OptimisationCriterion::MAXIMISE and rank < threshold)) {
                            if (p.severity() == RankingConstraintSeverity::PERMISSIVE) ++low_errors;
                            else if (p.severity() == RankingConstraintSeverity::CRITICAL) ++high_errors;
                        }
                    }
                    if (max_min_diff > 0) local_score = (rank-scalar_min_max[pw_idx].first)/max_min_diff;
                } else {
                    SizeType effective_dim = dim;
                    for (SizeType i=0; i<dim; ++i) {
                        auto max_min_diff = vector_min_max[pw_idx][i].second - vector_min_max[pw_idx][i].first;
                        auto rank = p.rank(input, entry.second.first, entry.second.second, i);
                        if (max_min_diff > 0) local_score = (rank - vector_min_max[pw_idx][i].first)/max_min_diff;
                        else --effective_dim;
                        if (effective_dim > 0) local_score/=effective_dim;
                    }
                }

                if (p.optimisation() == OptimisationCriterion::MAXIMISE) score += weight * local_score;
                else score -= weight*local_score;
            }
            result.push_back(TaskExecutionRanking(entry.first, score, low_errors, high_errors));
        }

        return result;
    }

  private:
    ParameterWeightsList<R> const _parameter_weights;
};

template<class R>
TaskRankingSpace<R> TaskRankingSpaceBuilder<R>::build() const {
    return TaskRankingSpace(_parameter_weights);
}

} // namespace Ariadne

#endif // ARIADNE_TASK_RANKING_SPACE_HPP
