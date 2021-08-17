#include "N_CORNER.h"
#include <cstdio>
#include <utility>
#include <vector>
#include <cassert>
#include <bitset>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iostream>

#define NUM_DIMS 2

/*
 *Code kindly provided by the authors from the following paper:

 *"Improving Spatial Data Processing by Clipping Minimum Bounding Boxes"
 * 2018 
 * Darius Sidlauskas; Sean Chester; Eleni Tzirita Zacharatou; Anastasia Ailamaki
*/

class Point_C
{
public:
    double coord[NUM_DIMS];

    inline bool operator<(const Point_C &rhs) const
    {
        for (uint32_t d = 0; d < NUM_DIMS; ++d)
        {
            if (coord[d] == rhs[d])
            {
                continue;
            }
            return coord[d] < rhs[d];
        }
        return false; // points are equal
    }

    inline bool operator==(const Point_C &rhs) const
    {
        for (uint32_t d = 0; d < NUM_DIMS; ++d)
            if (coord[d] != rhs[d])
                return false;

        return true;
    }

    inline bool operator!=(const Point_C &rhs) const
    {
        return !operator==(rhs);
    }

    inline double &operator[](const size_t idx)
    {
        assert(idx >= 0 && idx < NUM_DIMS);
        return coord[idx];
    };

    inline const double &operator[](const size_t idx) const
    {
        assert(idx >= 0 && idx < NUM_DIMS);
        return coord[idx];
    }
};

using PointPair = std::pair<Point_C, Point_C>; // for line segments.
using PointSet = std::vector<Point_C>;

static double euclidean_distance_squared(const Point_C &p, const Point_C &q)
{
    return std::inner_product(p.coord, p.coord + NUM_DIMS, q.coord, 0.0,
                              std::plus<double>(),
                              [](double pval, double qval)
                              {
                                  return (pval - qval) * (pval - qval);
                              });
}

static double euclidean_distance(const Point_C &p, const Point_C &q)
{
    return std::sqrt(euclidean_distance_squared(p, q));
}

static inline double dotProduct(const Point_C &p, const Point_C &q)
{
    return std::inner_product(p.coord, p.coord + NUM_DIMS, q.coord, 0.0);
}

static double vectorMagnitude(const Point_C &ref, const Point_C &p)
{
    return std::sqrt(
        std::inner_product(p.coord, p.coord + NUM_DIMS, ref.coord, 0.0,
                           std::plus<double>(),
                           [](double pval, double refval)
                           { return (pval - refval) * (pval - refval); }));
}

static double crossProduct(const Point_C &p, const Point_C &q)
{
    return p.coord[0] * q.coord[1] - q.coord[0] * p.coord[1];
}

static double areaOfPolygon(const PointSet &ps)
{
    return std::inner_product(ps.begin(), ps.end() - 1, ps.begin() + 1,
                              crossProduct(*(ps.end() - 1), *(ps.begin())),
                              std::plus<double>(),
                              crossProduct) /
           2;
}

static PointSet fiveCorners(const PointSet &convex_hull)
{
    PointSet result(5);
    const size_t num_points = convex_hull.size();

    if (num_points <= 5)
    {
        assert(num_points >= 3); // otherwise is it really a polygon?
        // then solution == input!
        return convex_hull;
    }

    // transform the input polygon from vertices into edges.
    std::vector<PointPair> edges(2 * num_points - 1);
    std::transform(convex_hull.begin(), convex_hull.end() - 1,
                   convex_hull.begin() + 1, edges.begin(),
                   [](const Point_C &p, const Point_C &q)
                   {
                       return PointPair(p, q);
                   });
    edges[num_points - 1] = PointPair(convex_hull[num_points - 1], convex_hull[0]);

    // duplicate the edge list for easier "wrapping around" (i.e., mod P)
    std::transform(convex_hull.begin(), convex_hull.end() - 1,
                   convex_hull.begin() + 1, edges.begin() + num_points,
                   [](const Point_C &p, const Point_C &q)
                   {
                       return PointPair(p, q);
                   });

    // finds the intersection point of line1 and line2
    auto intersect = [](const PointPair &line1, const PointPair &line2)
    {
        const double line1_rise = (line1.second[1] - line1.first[1]);
        const double line1_run = (line1.second[0] - line1.first[0]);
        const double line2_rise = (line2.second[1] - line2.first[1]);
        const double line2_run = (line2.second[0] - line2.first[0]);

        const double line1_y_intercept = line1.first[1] - line1_rise / line1_run * line1.first[0];
        const double line2_y_intercept = line2.first[1] - line2_rise / line2_run * line2.first[0];

        if ((line1_run == 0 && line2_run == 0) || (line1_rise / line1_run == line2_rise / line2_run))
        {
            // parallel lines. Return "infinite" intersection point.
            return Point_C{std::numeric_limits<double>::max(),
                           std::numeric_limits<double>::max()};
        }
        if (line1_run == 0) // vertical line causes div-by-zero
        {
            const double xval = line1.first[0];
            const double yval = line2_rise / line2_run * xval + line2_y_intercept;
            return Point_C{xval, yval};
        }
        else if (line2_run == 0)
        {
            const double xval = line2.first[0];
            const double yval = line1_rise / line1_run * xval + line1_y_intercept;
            return Point_C{xval, yval};
        }
        else
        {
            const double xval = (line2_y_intercept - line1_y_intercept) / (line1_rise / line1_run - line2_rise / line2_run);
            const double yval = line1_rise / line1_run * xval + line1_y_intercept;
            return Point_C{xval, yval};
        }
    };

    // calculates the area contained by the flush circumscribing polygon that
    // starts at i, ends at j, and contains edge mid.
    auto flush_added_area = [&intersect](const PointPair &i, const PointPair &j, const PointPair &mid,
                                         Point_C &new_i, Point_C &new_j)
    {
        const Point_C intersect_i_mid = intersect(i, mid);
        const Point_C intersect_mid_j = intersect(mid, j);

        if ( // lines are parallel
            intersect_i_mid[0] == std::numeric_limits<double>::max() || intersect_mid_j[0] == std::numeric_limits<double>::max()
            // || angle is too great thus intersection point is in wrong direction...
            || (euclidean_distance_squared(intersect_i_mid, i.first) < euclidean_distance_squared(intersect_i_mid, i.second)) || (euclidean_distance_squared(intersect_mid_j, j.second) < euclidean_distance_squared(intersect_mid_j, j.first)))
        {
            // could not close polygon. return impossibly high area
            // together with arbitrary point pair.
            return std::numeric_limits<double>::max();
        }
        else
        {
            new_i[0] = intersect_i_mid[0];
            new_i[1] = intersect_i_mid[1];
            new_j[0] = intersect_mid_j[0];
            new_j[1] = intersect_mid_j[1];

            return areaOfPolygon(
                PointSet{i.first, intersect_i_mid, intersect_mid_j, j.second});
        }
    };

    // calculates the area contained by the flush circumscribing polygon that
    // starts at i, ends at j, and contains two edges mid1 and mid2.
    auto twoflush_added_area = [&intersect](const PointPair &i, const PointPair &j, const PointPair &mid1,
                                            const PointPair &mid2, Point_C &new_i, Point_C &new_j, Point_C &new_mid)
    {
        const Point_C intersect_i_mid = intersect(i, mid1);
        const Point_C intersect_mids = intersect(mid1, mid2);
        const Point_C intersect_mid_j = intersect(mid2, j);

        if ( // lines are parallel
            intersect_i_mid[0] == std::numeric_limits<double>::max() || intersect_mid_j[0] == std::numeric_limits<double>::max() || intersect_mids[0] == std::numeric_limits<double>::max()
            // || angle is too great thus intersection point is in wrong direction...
            || (euclidean_distance_squared(intersect_i_mid, i.first) < euclidean_distance_squared(intersect_i_mid, i.second)) || (euclidean_distance_squared(intersect_mid_j, j.second) < euclidean_distance_squared(intersect_mid_j, j.first)) || (euclidean_distance_squared(intersect_mids, mid1.first) < euclidean_distance_squared(intersect_mids, mid1.second)))
        {
            // could not close polygon. return impossibly high area
            // together with arbitrary point pair.
            return std::numeric_limits<double>::max();
        }
        else
        {
            new_i[0] = intersect_i_mid[0];
            new_i[1] = intersect_i_mid[1];
            new_j[0] = intersect_mid_j[0];
            new_j[1] = intersect_mid_j[1];
            new_mid[0] = intersect_mids[0];
            new_mid[1] = intersect_mids[1];

            return areaOfPolygon(
                PointSet{i.first, intersect_i_mid, intersect_mids, intersect_mid_j, j.second});
        }
    };

    // calculates the area contained by the non-flush circumscribing polygon that
    // starts at i, ends at j, and passes through vertex mid.
    auto nonflush_added_area = [&intersect](const PointPair &i, const PointPair &j, const Point_C &mid,
                                            Point_C &new_i, Point_C &new_j)
    {
        const Point_C head = intersect(i, j);
        if (
            head[0] == std::numeric_limits<double>::max() || (euclidean_distance_squared(head, i.first) < euclidean_distance_squared(head, i.second)) || (euclidean_distance_squared(head, j.second) < euclidean_distance_squared(head, j.first)))
        {
            // does not close polygon. return impossibly high score
            // together with arbitrary point pair.
            return std::numeric_limits<double>::max();
        }

        const double mid_dist = euclidean_distance(head, mid);
        const Point_C mid_norm{(mid[0] - head[0]) / mid_dist,
                               (mid[1] - head[1]) / mid_dist};

        // double check that mid_dist is reasonable.
        // sometimes "near parallel" lines lead to divisions
        // that even "double" doesn't have sufficient precision for.
        if (mid_dist > 1000 * euclidean_distance(i.first, j.first))
        {
            return std::numeric_limits<double>::max();
        }

        // pick further point on line segment to avoid magnitudes close to 0.
        const Point_C &i_ref = (vectorMagnitude(i.first, head) < vectorMagnitude(i.second, head)
                                    ? i.second
                                    : i.first);
        const double i_mag = vectorMagnitude(i_ref, head);
        const Point_C i_norm{(i_ref[0] - head[0]) / i_mag,
                             (i_ref[1] - head[1]) / i_mag};

        const Point_C &j_ref = (vectorMagnitude(j.first, head) < vectorMagnitude(j.second, head)
                                    ? j.second
                                    : j.first);
        const double j_mag = vectorMagnitude(j_ref, head);
        const Point_C j_norm{(j_ref[0] - head[0]) / j_mag,
                             (j_ref[1] - head[1]) / j_mag};

        const double angle_i_mid = std::acos(dotProduct(i_norm, mid_norm));
        const double angle_j_mid = std::acos(dotProduct(mid_norm, j_norm));
        const double angle_i_j = std::acos(dotProduct(i_norm, j_norm));

        const double dist_new_j = mid_dist * 2 * sin(angle_i_mid) / sin(angle_i_j);
        const double dist_new_i = mid_dist * 2 * sin(angle_j_mid) / sin(angle_i_j);

        new_i[0] = head[0] + i_norm[0] * dist_new_i;
        new_i[1] = head[1] + i_norm[1] * dist_new_i;
        new_j[0] = head[0] + j_norm[0] * dist_new_j;
        new_j[1] = head[1] + j_norm[1] * dist_new_j;

        return areaOfPolygon(
            PointSet{i.first, new_i, new_j, j.second});
    };

    double global_best_area = std::numeric_limits<double>::max();
    using namespace std::placeholders;
    Point_C new_i, new_j, new_mid;
    for (size_t i = 0; i < num_points; ++i)
    {
        for (size_t j = i + 3; j < i + num_points - 2; ++j)
        {
            const double ij_quadrilateral_area = areaOfPolygon(
                PointSet{edges[i].first, edges[i].second, edges[j].first, edges[j].second});

            auto flush_ij = std::bind(
                twoflush_added_area, edges[i], edges[j], _1, _2,
                std::ref(new_i), std::ref(new_j), std::ref(new_mid));
            std::vector<std::pair<double, PointSet>> added_areas;
            for (auto mid1 = i + 1; mid1 < j - 1; ++mid1)
            {
                for (auto mid2 = mid1 + 1; mid2 < j; ++mid2)
                {
                    const double new_area = flush_ij(edges[mid1], edges[mid2]);
                    added_areas.push_back(std::make_pair(
                        new_area, PointSet{Point_C(new_i), Point_C(new_mid), Point_C(new_j)}));
                }
            }
            const auto min_flush = std::min_element(added_areas.begin(), added_areas.end());

            auto flush_ji = std::bind(
                flush_added_area, edges[j], edges[i], _1, std::ref(new_i), std::ref(new_j));
            auto nonflush_ji = std::bind(
                nonflush_added_area, edges[j], edges[i], _1, std::ref(new_i), std::ref(new_j));

            std::vector<std::pair<double, PointPair>> added_areas_bottom;
            for (size_t mid = j + 1; mid < i + num_points - 1; ++mid)
            {
                const double flush = flush_ji(edges[mid]);
                added_areas_bottom.push_back(std::make_pair(
                    flush, PointPair(Point_C(new_i), Point_C(new_j))));
                const double nonflush = nonflush_ji(edges[mid].second);
                added_areas_bottom.push_back(std::make_pair(
                    nonflush, PointPair(Point_C(new_i), Point_C(new_j))));
            }
            const double last_ji_flush = flush_ji(edges[i + num_points - 1]);
            added_areas_bottom.push_back(std::make_pair(
                last_ji_flush, PointPair(Point_C(new_i), Point_C(new_j))));
            const auto min_nonflush = std::min_element(
                added_areas_bottom.begin(), added_areas_bottom.end());

            const double total_area = min_flush->first + min_nonflush->first - ij_quadrilateral_area; // subtracting off the double-counted area between i/j.

            if (total_area < global_best_area)
            {
                global_best_area = total_area;
                result = PointSet{
                    min_flush->second[0], min_flush->second[1], min_flush->second[2],
                    min_nonflush->second.first, min_nonflush->second.second};
            }
        }
    }

    return result;
}

static PointSet fourCorners(const PointSet &convex_hull)
{
    PointSet result(4);
    const size_t num_points = convex_hull.size();

    if (num_points <= 4)
    {
        assert(num_points >= 3); // otherwise is it really a polygon?
        // then solution == input!
        return convex_hull;
    }

    // transform the input polygon from vertices into edges.
    std::vector<PointPair> edges(2 * num_points - 1);
    std::transform(convex_hull.begin(), convex_hull.end() - 1,
                   convex_hull.begin() + 1, edges.begin(),
                   [](const Point_C &p, const Point_C &q)
                   {
                       return PointPair(p, q);
                   });
    edges[num_points - 1] = PointPair(convex_hull[num_points - 1], convex_hull[0]);

    // duplicate the edge list for easier "wrapping around" (i.e., mod P)
    std::transform(convex_hull.begin(), convex_hull.end() - 1,
                   convex_hull.begin() + 1, edges.begin() + num_points,
                   [](const Point_C &p, const Point_C &q)
                   {
                       return PointPair(p, q);
                   });

    // finds the intersection point of line1 and line2
    auto intersect = [](const PointPair &line1, const PointPair &line2)
    {
        const double line1_rise = (line1.second[1] - line1.first[1]);
        const double line1_run = (line1.second[0] - line1.first[0]);
        const double line2_rise = (line2.second[1] - line2.first[1]);
        const double line2_run = (line2.second[0] - line2.first[0]);

        const double line1_y_intercept = line1.first[1] - line1_rise / line1_run * line1.first[0];
        const double line2_y_intercept = line2.first[1] - line2_rise / line2_run * line2.first[0];

        if ((line1_run == 0 && line2_run == 0) || (line1_rise / line1_run == line2_rise / line2_run))
        {
            // parallel lines. Return "infinite" intersection point.
            return Point_C{std::numeric_limits<double>::max(),
                           std::numeric_limits<double>::max()};
        }
        if (line1_run == 0) // vertical line causes div-by-zero
        {
            const double xval = line1.first[0];
            const double yval = line2_rise / line2_run * xval + line2_y_intercept;
            return Point_C{xval, yval};
        }
        else if (line2_run == 0)
        {
            const double xval = line2.first[0];
            const double yval = line1_rise / line1_run * xval + line1_y_intercept;
            return Point_C{xval, yval};
        }
        else
        {
            const double xval = (line2_y_intercept - line1_y_intercept) / (line1_rise / line1_run - line2_rise / line2_run);
            const double yval = line1_rise / line1_run * xval + line1_y_intercept;
            return Point_C{xval, yval};
        }
    };

    // calculates the area contained by the flush circumscribing polygon that
    // starts at i, ends at j, and contains edge mid.
    auto flush_added_area = [&intersect](const PointPair &i, const PointPair &j, const PointPair &mid,
                                         Point_C &new_i, Point_C &new_j)
    {
        const Point_C intersect_i_mid = intersect(i, mid);
        const Point_C intersect_mid_j = intersect(mid, j);

        if ( // lines are parallel
            intersect_i_mid[0] == std::numeric_limits<double>::max() || intersect_mid_j[0] == std::numeric_limits<double>::max()
            // || angle is too great thus intersection point is in wrong direction...
            || (euclidean_distance_squared(intersect_i_mid, i.first) < euclidean_distance_squared(intersect_i_mid, i.second)) || (euclidean_distance_squared(intersect_mid_j, j.second) < euclidean_distance_squared(intersect_mid_j, j.first)))
        {
            // could not close polygon. return impossibly high area
            // together with arbitrary point pair.
            return std::numeric_limits<double>::max();
        }
        else
        {
            new_i[0] = intersect_i_mid[0];
            new_i[1] = intersect_i_mid[1];
            new_j[0] = intersect_mid_j[0];
            new_j[1] = intersect_mid_j[1];

            return areaOfPolygon(
                PointSet{i.first, intersect_i_mid, intersect_mid_j, j.second});
        }
    };

    // calculates the area contained by the non-flush circumscribing polygon that
    // starts at i, ends at j, and passes through vertex mid.
    auto nonflush_added_area = [&intersect](const PointPair &i, const PointPair &j, const Point_C &mid,
                                            Point_C &new_i, Point_C &new_j)
    {
        const Point_C head = intersect(i, j);
        if (
            head[0] == std::numeric_limits<double>::max() || (euclidean_distance_squared(head, i.first) < euclidean_distance_squared(head, i.second)) || (euclidean_distance_squared(head, j.second) < euclidean_distance_squared(head, j.first)))
        {
            // does not close polygon. return impossibly high score
            // together with arbitrary point pair.
            return std::numeric_limits<double>::max();
        }

        const double mid_dist = euclidean_distance(head, mid);
        const Point_C mid_norm{(mid[0] - head[0]) / mid_dist,
                               (mid[1] - head[1]) / mid_dist};

        // double check that mid_dist is reasonable.
        // sometimes "near parallel" lines lead to divisions
        // that even "double" doesn't have sufficient precision for.
        if (mid_dist > 1000 * euclidean_distance(i.first, j.first))
        {
            return std::numeric_limits<double>::max();
        }

        // pick further point on line segment to avoid magnitudes close to 0.
        const Point_C &i_ref = (vectorMagnitude(i.first, head) < vectorMagnitude(i.second, head)
                                    ? i.second
                                    : i.first);
        const double i_mag = vectorMagnitude(i_ref, head);
        const Point_C i_norm{(i_ref[0] - head[0]) / i_mag,
                             (i_ref[1] - head[1]) / i_mag};

        const Point_C &j_ref = (vectorMagnitude(j.first, head) < vectorMagnitude(j.second, head)
                                    ? j.second
                                    : j.first);
        const double j_mag = vectorMagnitude(j_ref, head);
        const Point_C j_norm{(j_ref[0] - head[0]) / j_mag,
                             (j_ref[1] - head[1]) / j_mag};

        const double angle_i_mid = std::acos(dotProduct(i_norm, mid_norm));
        const double angle_j_mid = std::acos(dotProduct(mid_norm, j_norm));
        const double angle_i_j = std::acos(dotProduct(i_norm, j_norm));

        const double dist_new_j = mid_dist * 2 * sin(angle_i_mid) / sin(angle_i_j);
        const double dist_new_i = mid_dist * 2 * sin(angle_j_mid) / sin(angle_i_j);

        new_i[0] = head[0] + i_norm[0] * dist_new_i;
        new_i[1] = head[1] + i_norm[1] * dist_new_i;
        new_j[0] = head[0] + j_norm[0] * dist_new_j;
        new_j[1] = head[1] + j_norm[1] * dist_new_j;

        return areaOfPolygon(
            PointSet{i.first, new_i, new_j, j.second});
    };

    double global_best_area = std::numeric_limits<double>::max();
    using namespace std::placeholders;
    Point_C new_i, new_j;
    for (size_t i = 0; i < num_points; ++i)
    {
        for (size_t j = i + 2; j < i + num_points - 2; ++j)
        {
            const double ij_quadrilateral_area = areaOfPolygon(
                PointSet{edges[i].first, edges[i].second, edges[j].first, edges[j].second});

            auto flush_ij = std::bind(
                flush_added_area, edges[i], edges[j], _1, std::ref(new_i), std::ref(new_j));
            std::vector<std::pair<double, PointPair>> added_areas;
            for (auto mid = i + 1; mid < j; ++mid)
            {
                const double new_area = flush_ij(edges[mid]);
                added_areas.push_back(std::make_pair(
                    new_area, PointPair(Point_C(new_i), Point_C(new_j))));
            }
            const auto min_flush = std::min_element(added_areas.begin(), added_areas.end());

            auto flush_ji = std::bind(
                flush_added_area, edges[j], edges[i], _1, std::ref(new_i), std::ref(new_j));
            auto nonflush_ji = std::bind(
                nonflush_added_area, edges[j], edges[i], _1, std::ref(new_i), std::ref(new_j));

            std::vector<std::pair<double, PointPair>> added_areas_bottom;
            for (size_t mid = j + 1; mid < i + num_points - 1; ++mid)
            {
                const double flush = flush_ji(edges[mid]);
                added_areas_bottom.push_back(std::make_pair(
                    flush, PointPair(Point_C(new_i), Point_C(new_j))));
                const double nonflush = nonflush_ji(edges[mid].second);
                added_areas_bottom.push_back(std::make_pair(
                    nonflush, PointPair(Point_C(new_i), Point_C(new_j))));
            }
            const double last_ji_flush = flush_ji(edges[i + num_points - 1]);
            added_areas_bottom.push_back(std::make_pair(
                last_ji_flush, PointPair(Point_C(new_i), Point_C(new_j))));
            const auto min_nonflush = std::min_element(
                added_areas_bottom.begin(), added_areas_bottom.end());

            const double total_area = min_flush->first + min_nonflush->first - ij_quadrilateral_area; // subtracting off the double-counted area between i/j.

            if (total_area < global_best_area)
            {
                global_best_area = total_area;
                result = PointSet{min_flush->second.first, min_flush->second.second,
                                  min_nonflush->second.first, min_nonflush->second.second};
            }
        }
    }

    return result;
}

static PointSet vector_to_point_set(int n, ApproxPoint *P)
{
    PointSet set;
    for (int i = 0; i < n; ++i)
    {
        Point_C t;
        t.coord[0] = P[i].x;
        t.coord[1] = P[i].y;
        set.push_back(t);
    }
    return set;
}

extern "C"
{

    N_CORNER five_corner_calculate(int n, ApproxPoint *P)
    {
        PointSet _5corner = fiveCorners(vector_to_point_set(n, P));
        ApproxPoint *vec = new ApproxPoint[_5corner.size()];
        PointSet::iterator it = _5corner.begin();
        int i = 0;
        while (it != _5corner.end())
        {
            ApproxPoint p;
            p.x = (*it).coord[0];
            p.y = (*it).coord[1];
            vec[i] = p;
            ++i;
            ++it;
        }
        N_CORNER n_corner;
        n_corner.p = vec;
        n_corner.num = _5corner.size();
        return n_corner;
    }

    N_CORNER four_corner_calculate(int n, ApproxPoint *P)
    {
        PointSet _4corner = fourCorners(vector_to_point_set(n, P));
        ApproxPoint *vec = new ApproxPoint[_4corner.size()];
        PointSet::iterator it = _4corner.begin();
        int i = 0;
        while (it != _4corner.end())
        {
            ApproxPoint p;
            p.x = (*it).coord[0];
            p.y = (*it).coord[1];
            vec[i] = p;
            ++i;
            ++it;
        }
        N_CORNER n_corner;
        n_corner.p = vec;
        n_corner.num = _4corner.size();
        return n_corner;
    }

    ApproxPoint* n_corner_allocate_points(int n)
    {
        return new ApproxPoint[n];
    }

    void n_corner_free(N_CORNER *n_corner)
    {
        delete[] n_corner->p;
    }
}
