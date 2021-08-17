#include "RMBR.h"
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Polygon_2<Kernel> Polygon_2;

extern "C" {

    RMBR rmbr_create(ApproxPoint a, ApproxPoint b, ApproxPoint c, ApproxPoint d) {
        Polygon_2* polygon = new Polygon_2();
        polygon->push_back(Point_2(a.x, a.y));
        polygon->push_back(Point_2(b.x, b.y));
        polygon->push_back(Point_2(c.x, c.y));
        polygon->push_back(Point_2(d.x, d.y));
        return (RMBR) polygon;
    }

    void rmbr_free(RMBR rmbr) {
        delete ((Polygon_2*) rmbr);
    }

    bool rmbr_is_in_exterior(RMBR rmbr, ApproxPoint p) {
        Polygon_2* poly = ((Polygon_2*) rmbr);
        return poly->has_on_unbounded_side(Point_2(p.x, p.y));
    }

    bool rmbr_is_in_boundary(RMBR rmbr, ApproxPoint p) {
        Polygon_2* poly = ((Polygon_2*) rmbr);
        return poly->has_on_boundary(Point_2(p.x, p.y));
    }

    bool rmbr_is_in_interior(RMBR rmbr, ApproxPoint p) {
        Polygon_2* poly = ((Polygon_2*) rmbr);
        return poly->has_on_bounded_side(Point_2(p.x, p.y));
    }

    int rmbr_region(RMBR rmbr, ApproxPoint p) {
        Polygon_2* poly = ((Polygon_2*) rmbr);
        return poly->bounded_side(Point_2(p.x, p.y));
    }

    RMBR rmbr_calculate(int n, ApproxPoint* P) {
        Polygon_2 poly;
        Polygon_2* rect = new Polygon_2();
        for (int i = 0; i < n; i++) {
            poly.push_back(Point_2((P[i]).x, (P[i]).y));
        }
        CGAL::min_rectangle_2(poly.vertices_begin(), poly.vertices_end(), std::back_inserter(*rect));
        return (RMBR) rect;
    }

    ApproxPoint *rmbr_get_points(RMBR rmbr) {
        Polygon_2* poly = ((Polygon_2*) rmbr);
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