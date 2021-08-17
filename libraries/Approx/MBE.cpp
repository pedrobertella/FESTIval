#include "MBE.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Min_ellipse_2/Optimisation_ellipse_2.h>
#include <iostream>

typedef CGAL::Exact_rational NT;
typedef CGAL::Cartesian<NT> K;
typedef CGAL::Point_2<K> Point_2;
typedef CGAL::Min_ellipse_2_traits_2<K> Traits;
typedef CGAL::Min_ellipse_2<Traits> Min_ellipse;

extern "C" {

    MBE mbe_create(int n, ApproxPoint *P) {
        CGAL::Optimisation_ellipse_2<K>* ellipse = new CGAL::Optimisation_ellipse_2<K>();
        if (n == 0) {
            ellipse->set();
        }
        if (n >= 1) {
            ellipse->set(Point_2((P[0]).x, (P[0]).y));
        }
        if (n >= 2) {
            ellipse->set(Point_2((P[0]).x, (P[0]).y), Point_2((P[1]).x, (P[1]).y));
        }
        if (n >= 3) {
            ellipse->set(Point_2((P[0]).x, (P[0]).y), Point_2((P[1]).x, (P[1]).y), Point_2((P[2]).x, (P[2]).y));
        }
        if (n >= 4) {
            ellipse->set(Point_2((P[0]).x, (P[0]).y), Point_2((P[1]).x, (P[1]).y), Point_2((P[2]).x, (P[2]).y), Point_2((P[3]).x, (P[3]).y));
        }
        if (n == 5) {
            ellipse->set(Point_2((P[0]).x, (P[0]).y), Point_2((P[1]).x, (P[1]).y), Point_2((P[2]).x, (P[2]).y), Point_2((P[3]).x, (P[3]).y), Point_2((P[4]).x, (P[4]).y));
        }
        return (MBE) ellipse;
    }

    void mbe_free(MBE mbe) {
        delete ((CGAL::Optimisation_ellipse_2<K>*)mbe);
    }

    bool mbe_is_in_exterior(MBE mbe, ApproxPoint p) {
        CGAL::Optimisation_ellipse_2<K>* ellipse = ((CGAL::Optimisation_ellipse_2<K>*)mbe);
        return ellipse->has_on_unbounded_side(Point_2(p.x, p.y));
    }

    bool mbe_is_in_boundary(MBE mbe, ApproxPoint p) {
        CGAL::Optimisation_ellipse_2<K>* ellipse = ((CGAL::Optimisation_ellipse_2<K>*)mbe);
        return ellipse->has_on_boundary(Point_2(p.x, p.y));
    }

    bool mbe_is_in_interior(MBE mbe, ApproxPoint p) {
        CGAL::Optimisation_ellipse_2<K>* ellipse = ((CGAL::Optimisation_ellipse_2<K>*)mbe);
        return ellipse->has_on_bounded_side(Point_2(p.x, p.y));
    }

    int mbe_region(MBE mbe, ApproxPoint p) {
        CGAL::Optimisation_ellipse_2<K>* ellipse = ((CGAL::Optimisation_ellipse_2<K>*)mbe);
        return ellipse->bounded_side(Point_2(p.x, p.y));
    }

    MBE mbe_calculate(int n, ApproxPoint* P) {
        Point_2 *points = new Point_2[n];
        for (int i = 0; i < n; i++) {
            points[i] = Point_2((P[i]).x, (P[i]).y);
        }
        Min_ellipse me(points, points + n, true);
        CGAL::Optimisation_ellipse_2<K>* ellipse = new CGAL::Optimisation_ellipse_2<K>();
        *ellipse = me.ellipse();
        delete[] points;
        return (MBE) ellipse;
    }

    ApproxPoint *mbe_get_points(int *n, MBE mbe) {
        CGAL::Optimisation_ellipse_2<K>* ellipse = ((CGAL::Optimisation_ellipse_2<K>*)mbe);
        *n = ellipse->number_of_boundary_points();
        ApproxPoint *vec = new ApproxPoint[*n];
        if(*n>=1){
            vec[0].x = CGAL::to_double(ellipse->boundary_point1.x());
            vec[0].y = CGAL::to_double(ellipse->boundary_point1.y());
        }
        if(*n>=2){
            vec[1].x = CGAL::to_double(ellipse->boundary_point2.x());
            vec[1].y = CGAL::to_double(ellipse->boundary_point2.y());
        }
        if(*n>=3){
            vec[2].x = CGAL::to_double(ellipse->boundary_point3.x());
            vec[2].y = CGAL::to_double(ellipse->boundary_point3.y());
        }
        if(*n>=4){
            vec[3].x = CGAL::to_double(ellipse->boundary_point4.x());
            vec[3].y = CGAL::to_double(ellipse->boundary_point4.y());
        }
        if(*n==5){
            vec[4].x = CGAL::to_double(ellipse->boundary_point5.x());
            vec[4].y = CGAL::to_double(ellipse->boundary_point5.y());
        }
        return vec;
    }
}