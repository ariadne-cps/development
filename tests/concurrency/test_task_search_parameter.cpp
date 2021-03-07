/***************************************************************************
 *            test_task_search_parameter.cpp
 *
 *  Copyright  2008-20  Luca Geretti
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

#include "symbolic/expression_set.hpp"
#include "configuration/configuration_search_point.hpp"
#include "configuration/configuration_search_space.hpp"
#include "concurrency/task_execution_ranking.hpp"

#include "../test.hpp"

using namespace Ariadne;

class TestTaskSearchParameter {
  public:

    Void test_task_parameter_creation() {
        ConfigurationSearchParameter p(ConfigurationPropertyPath("use_subdivisions"), false, List<int>({0, 1}));
        ARIADNE_TEST_PRINT(p);
    }

    Void test_task_parameter_randomise() {
        ConfigurationSearchParameter p(ConfigurationPropertyPath("use_subdivisions"), false, List<int>({0, 1}));
        ARIADNE_TEST_PRINT(p);
        List<int> values;
        for (Nat i=0; i<16; ++i) values.push_back(p.random_value());
        ARIADNE_TEST_PRINT(values);
    }

    Void test_metric_task_parameter_shift() {
        ConfigurationSearchParameter metric(ConfigurationPropertyPath("sweep_threshold"), true, List<int>({8, 9, 10, 11}));
        ARIADNE_TEST_EQUALS(metric.shifted_value_from(8),9);
        ARIADNE_TEST_EQUALS(metric.shifted_value_from(11),10);
        auto from_1 = metric.shifted_value_from(10);
        ARIADNE_TEST_ASSERT(from_1 == 9 or from_1 == 11);
    }

    Void test_parameter_space() {
        ConfigurationPropertyPath use_subdivisions("use_subdivisions");
        ConfigurationPropertyPath sweep_threshold("sweep_threshold");
        ConfigurationSearchParameter bp(use_subdivisions, false, List<int>({0, 1}));
        ConfigurationSearchParameter mp(sweep_threshold, true, List<int>({3, 4, 5}));
        ConfigurationSearchSpace space({bp, mp});
        ARIADNE_TEST_PRINT(space);
        ARIADNE_TEST_PRINT(space.parameters());
        ARIADNE_TEST_EQUALS(space.dimension(),2);
        ARIADNE_TEST_EQUALS(space.index(bp),1);
        ARIADNE_TEST_EQUALS(space.index(mp),0);
        ARIADNE_TEST_EQUALS(space.total_points(),6);
    }

    Void test_parameter_point_creation() {
        ConfigurationPropertyPath use_subdivisions("use_subdivisions");
        ConfigurationPropertyPath sweep_threshold("sweep_threshold");
        ConfigurationSearchParameter bp(use_subdivisions, false, List<int>({0, 1}));
        ConfigurationSearchParameter mp(sweep_threshold, true, List<int>({3, 4, 5}));
        ConfigurationSearchSpace space({bp, mp});
        ARIADNE_TEST_PRINT(space.initial_point());
        Map<ConfigurationPropertyPath,int> bindings = {{use_subdivisions,1},{sweep_threshold,5}};
        ARIADNE_TEST_PRINT(bindings);
        ConfigurationSearchPoint point = space.make_point(bindings);
        ARIADNE_TEST_PRINT(point);
        ARIADNE_TEST_PRINT(point.space());
    }

    Void test_parameter_point_equality() {
        ConfigurationPropertyPath use_subdivisions("use_subdivisions");
        ConfigurationPropertyPath sweep_threshold("sweep_threshold");
        ConfigurationSearchParameter bp(use_subdivisions, false, List<int>({0, 1}));
        ConfigurationSearchParameter mp(sweep_threshold, true, List<int>({3, 4, 5}));
        ConfigurationSearchSpace space({bp, mp});

        ConfigurationSearchPoint point1 = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 2}});
        ConfigurationSearchPoint point2 = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 2}});
        ConfigurationSearchPoint point3 = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 3}});
        ARIADNE_TEST_EQUAL(point1,point2);
        ARIADNE_TEST_NOT_EQUAL(point1,point3);
    }

    Void test_parameter_point_distance() {
        ConfigurationPropertyPath use_subdivisions("use_subdivisions");
        ConfigurationPropertyPath sweep_threshold("sweep_threshold");
        ConfigurationSearchParameter bp(use_subdivisions, false, List<int>({0, 1}));
        ConfigurationSearchParameter mp(sweep_threshold, true, List<int>({3, 4, 5}));
        ConfigurationSearchSpace space({bp, mp});

        ConfigurationSearchPoint point1 = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 2}});
        ConfigurationSearchPoint point2 = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 2}});
        ConfigurationSearchPoint point3 = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 3}});
        ConfigurationSearchPoint point4 = space.make_point({{use_subdivisions, 0}, {sweep_threshold, 5}});
        ARIADNE_TEST_EQUALS(point1.distance(point2),0);
        ARIADNE_TEST_EQUALS(point1.distance(point3),1);
        ARIADNE_TEST_EQUALS(point3.distance(point4),3);
    }

    Void test_parameter_point_adjacent_shift() {
        ConfigurationPropertyPath use_subdivisions("use_subdivisions");
        ConfigurationPropertyPath sweep_threshold("sweep_threshold");
        ConfigurationSearchParameter bp(use_subdivisions, false, List<int>({0, 1}));
        ConfigurationSearchParameter mp(sweep_threshold, true, List<int>({3, 4, 5}));
        ConfigurationSearchSpace space({bp, mp});

        ConfigurationSearchPoint point1 = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 5}});
        ARIADNE_TEST_PRINT(point1);
        auto point2 = point1.make_adjacent_shifted();
        ARIADNE_TEST_PRINT(point2);
        ARIADNE_TEST_NOT_EQUAL(point1,point2);
    }

    Void test_parameter_point_random_shift() {
        ConfigurationPropertyPath use_subdivisions("use_subdivisions");
        ConfigurationPropertyPath sweep_threshold("sweep_threshold");
        ConfigurationSearchParameter bp(use_subdivisions, false, List<int>({0, 1}));
        ConfigurationSearchParameter mp(sweep_threshold, true, List<int>({3, 4, 5, 6, 7}));
        ConfigurationSearchSpace space({bp, mp});

        ConfigurationSearchPoint point = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 5}});
        auto points = point.make_random_shifted(1);
        ARIADNE_TEST_EQUALS(points.size(),1);
        points = point.make_random_shifted(3);
        ARIADNE_TEST_EQUALS(points.size(),3);
        points = point.make_random_shifted(space.total_points());
        ARIADNE_TEST_EQUALS(points.size(),space.total_points());
    }

    Void test_parameter_point_adjacent_set_shift() {
        ConfigurationPropertyPath use_subdivisions("use_subdivisions");
        ConfigurationPropertyPath sweep_threshold("sweep_threshold");
        ConfigurationSearchParameter bp(use_subdivisions, false, List<int>({0, 1}));
        ConfigurationSearchParameter mp(sweep_threshold, true, List<int>({3, 4, 5, 6, 7, 8}));
        ConfigurationSearchSpace space({bp, mp});
        ARIADNE_TEST_PRINT(space.total_points());

        ConfigurationSearchPoint point = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 5}});
        Set<ConfigurationSearchPoint> points = point.make_random_shifted(3);
        ARIADNE_TEST_PRINT(points);
        auto all_points = make_extended_set_by_shifting(points, 5);
        ARIADNE_TEST_EQUALS(all_points.size(),5);
        ARIADNE_TEST_PRINT(all_points);

        ConfigurationSearchPoint point1 = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 8}});
        ConfigurationSearchPoint point2 = space.make_point({{use_subdivisions, 0}, {sweep_threshold, 3}});
        Set<ConfigurationSearchPoint> border_points = {point1, point2};
        ARIADNE_TEST_PRINT(border_points);
        ARIADNE_PRINT_TEST_COMMENT("Checking maximum number of single shift points including the original border points");
        auto six_points = make_extended_set_by_shifting(border_points, 6);
        ARIADNE_TEST_PRINT(six_points);
        ARIADNE_PRINT_TEST_COMMENT("Checking 1 point over the number of possible adjacent shiftings");
        auto seven_points = make_extended_set_by_shifting(border_points, 7);
        ARIADNE_TEST_PRINT(seven_points);
        ARIADNE_PRINT_TEST_COMMENT("Checking up to the maximum number");
        auto twelve_points = make_extended_set_by_shifting(border_points, 12);
        ARIADNE_TEST_PRINT(twelve_points);
        ARIADNE_PRINT_TEST_COMMENT("Checking 1 point over the maximum number");
        ARIADNE_TEST_FAIL(make_extended_set_by_shifting(points, point.space().total_points()+1));
    }

    Void test_parameter_point_ranking() {
        ConfigurationPropertyPath use_subdivisions("use_subdivisions");
        ConfigurationPropertyPath sweep_threshold("sweep_threshold");
        ConfigurationSearchParameter bp(use_subdivisions, false, List<int>({0, 1}));
        ConfigurationSearchParameter mp(sweep_threshold, true, List<int>({3, 4, 5, 6, 7}));
        ConfigurationSearchSpace space({bp, mp});

        ConfigurationSearchPoint point1 = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 2}});
        ConfigurationSearchPoint point2 = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 2}});
        ConfigurationSearchPoint point3 = space.make_point({{use_subdivisions, 1}, {sweep_threshold, 3}});
        ConfigurationSearchPoint point4 = space.make_point({{use_subdivisions, 0}, {sweep_threshold, 4}});
        TaskExecutionRanking a1(point1, 2, 0, 0);
        TaskExecutionRanking a2(point2, 4, 1, 0);
        TaskExecutionRanking a3(point3, 4, 0, 0);
        TaskExecutionRanking a4(point4, 4, 0, 1);
        Set<TaskExecutionRanking> as = {a1, a2, a3, a4};

        ARIADNE_TEST_PRINT(as);
        ARIADNE_TEST_ASSERT(a2 < a1);
        ARIADNE_TEST_ASSERT(a1 < a3);
        ARIADNE_TEST_ASSERT(a2 < a3);
        ARIADNE_TEST_ASSERT(a4 < a2);
        ARIADNE_TEST_ASSERT(a4 < a3);
    }

    Void test() {
        ARIADNE_TEST_CALL(test_task_parameter_creation());
        ARIADNE_TEST_CALL(test_task_parameter_randomise());
        ARIADNE_TEST_CALL(test_metric_task_parameter_shift());
        ARIADNE_TEST_CALL(test_parameter_space());
        ARIADNE_TEST_CALL(test_parameter_point_creation());
        ARIADNE_TEST_CALL(test_parameter_point_equality());
        ARIADNE_TEST_CALL(test_parameter_point_distance());
        ARIADNE_TEST_CALL(test_parameter_point_adjacent_shift());
        ARIADNE_TEST_CALL(test_parameter_point_random_shift());
        ARIADNE_TEST_CALL(test_parameter_point_adjacent_set_shift());
        ARIADNE_TEST_CALL(test_parameter_point_ranking());
    }
};

int main() {
    TestTaskSearchParameter().test();
    return ARIADNE_TEST_FAILURES;
}
