#include "RMBP.h"
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Polygon_2<Kernel> Polygon_2;

extern "C" {

    RMBP rmbp_create(ApproxPoint a, ApproxPoint b, ApproxPoint c, ApproxPoint d) {
        Polygon_2* polygon = new Polygon_2();
        polygon->push_back(Point_2(a.x, a.y));
        polygon->push_back(Point_2(b.x, b.y));
        polygon->push_back(Point_2(c.x, c.y));
        polygon->push_back(Point_2(d.x, d.y));
        return (RMBP) polygon;
    }

    void rmbp_free(RMBP rmbp) {
        delete ((Polygon_2*) rmbp);
    }

    bool rmbp_is_in_exterior(RMBP rmbp, ApproxPoint p) {
        Polygon_2* poly = ((Polygon_2*) rmbp);
        return poly->has_on_unbounded_side(Point_2(p.x, p.y));
    }

    bool rmbp_is_in_boundary(RMBP rmbp, ApproxPoint p) {
        Polygon_2* poly = ((Polygon_2*) rmbp);
        return poly->has_on_boundary(Point_2(p.x, p.y));
    }

    bool rmbp_is_in_interior(RMBP rmbp, ApproxPoint p) {
        Polygon_2* poly = ((Polygon_2*) rmbp);
        return poly->has_on_bounded_side(Point_2(p.x, p.y));
    }

    int rmbp_region(RMBP rmbp, ApproxPoint p) {
        Polygon_2* poly = ((Polygon_2*) rmbp);
        return poly->bounded_side(Point_2(p.x, p.y));
    }

    RMBP rmbp_calculate(int n, ApproxPoint* P) {
        Polygon_2 poly;
        Polygon_2* rect = new Polygon_2();
        for (int i = 0; i < n; i++) {
            poly.push_back(Point_2((P[i]).x, (P[i]).y));
        }
        CGAL::min_parallelogram_2(poly.vertices_begin(), poly.vertices_end(), std::back_inserter(*rect));
        return (RMBP) rect;
    }

    ApproxPoint *rmbp_get_points(RMBP rmbp) {
        Polygon_2* poly = ((Polygon_2*) rmbp);
        ApproxPoint *vec = new ApproxPoint[4];
        int count = 0;
        Polygon_2::Vertex_const_iterator i = poly->vertices_begin();
        while (i != poly->vertices_end()) {
            vec[count].x = (*i).x();
            vec[count].y = (*i).y();
            ++i;
            ++count;
        }
        return vec;
    }
}