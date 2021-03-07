/***************************************************************************
 *            concurrency/concurrency_manager.cpp
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

#include "../utility/macros.hpp"
#include "../concurrency/concurrency_manager.hpp"

namespace Ariadne {

ConcurrencyManager::ConcurrencyManager() : _maximum_concurrency(std::thread::hardware_concurrency()), _concurrency(1) {}

SizeType ConcurrencyManager::maximum_concurrency() const {
    return _maximum_concurrency;
}

SizeType ConcurrencyManager::concurrency() const {
    return _concurrency;
}

void ConcurrencyManager::set_concurrency(SizeType value) {
    ARIADNE_PRECONDITION(value <= _maximum_concurrency and value > 0);
    std::lock_guard<std::mutex> lock(_data_mutex);
    _concurrency = value;
}

List<TaskExecutionRanking> ConcurrencyManager::best_rankings() const {
    return _best_rankings;
}

void ConcurrencyManager::append_best_ranking(TaskExecutionRanking const& ranking) {
    std::lock_guard<std::mutex> lock(_data_mutex);
    _best_rankings.push_back(ranking);
}

void ConcurrencyManager::clear_best_rankings() {
    _best_rankings.clear();
}

List<int> ConcurrencyManager::optimal_point() const {
    List<int> result;
    if (not _best_rankings.empty()) {
        auto space = _best_rankings.front().point().space();
        auto dimension = space.dimension();

        List<Map<int,SizeType>> frequencies;
        for (SizeType i=0; i<dimension; ++i) frequencies.push_back(Map<int,SizeType>());
        for (auto ranking : _best_rankings) {
            auto coordinates = ranking.point().coordinates();
            for (SizeType i=0; i<dimension; ++i) {
                auto iter = frequencies[i].find(coordinates[i]);
                if (iter == frequencies[i].end()) frequencies[i].insert(make_pair(coordinates[i],1));
                else frequencies[i].at(coordinates[i])++;
            }
        }

        for (SizeType i=0; i<dimension; ++i) {
            auto freq_it = frequencies[i].begin();
            int best_value = freq_it->first;
            SizeType best_frequency = freq_it->second;
            ++freq_it;
            while(freq_it != frequencies[i].end()) {
                if (freq_it->second >best_frequency) {
                    best_value = freq_it->first;
                    best_frequency = freq_it->second;
                }
                ++freq_it;
            }
            result.push_back(best_value);
        }
    }
    return result;
}

void ConcurrencyManager::print_best_rankings() const {
    if (not _best_rankings.empty()) {
        std::ofstream file;
        file.open("points.m");
        auto space = _best_rankings.front().point().space();
        auto dimension = space.dimension();
        auto size = _best_rankings.size();
        file << "x = [1:" << size << "];\n";
        Map<SizeType,List<int>> values;
        List<int> soft_failures;
        for (SizeType i=0; i<dimension; ++i) values.insert(make_pair(i,List<int>()));
        for (auto ranking : _best_rankings) {
            auto point = ranking.point();
            for (SizeType i=0; i<dimension; ++i) values.at(i).push_back(point.coordinates()[i]);
            soft_failures.push_back(ranking.permissive_failures());
        }
        file << "figure(1);\n";
        file << "hold on;\n";
        for (SizeType i=0; i<dimension; ++i) {
            file << "y" << i << " = " << values.at(i) << ";\n";
            auto name = space.parameters()[i].path().last();
            std::replace(name.begin(), name.end(), '_', ' ');
            file << "plot(x,y" << i << ",'DisplayName','" << name << "');\n";
        }
        file << "yf = " << soft_failures << ";\n";
        file << "plot(x,yf,'DisplayName','(soft failures)');\n";

        file << "legend;\n";
        file << "hold off;\n";
        file.close();
    }
}

}