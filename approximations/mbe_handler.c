#include <string.h>
#include <liblwgeom.h>
#include "mbe_handler.h"
#include "../main/lwgeom_handler.h"
#include "../main/math_util.h"
#include "unistd.h"
#include "../main/log_messages.h"

static bool intersect(const MBE mbe, const BBox *bbox);
static bool overlap(const MBE mbe, const BBox *bbox);
static bool meet(const MBE mbe, const BBox *bbox);
static bool coveredBy(const MBE mbe, const BBox *bbox);
static bool inside(const MBE mbe, const BBox *bbox);
static bool covers(const MBE mbe, const BBox *bbox);
static bool contains(const MBE mbe, const BBox *bbox);
static bool inside_or_coveredBy(const MBE mbe, const BBox *bbox);
static bool contain_or_covers(const MBE mbe, const BBox *bbox);
static bool equal(const MBE mbe, const BBox *bbox);
static double **discretize_line(int *n, const double *p1, const double *p2, double step, int dim);
static void free_point_matrix(double **m, int n);
static bool region_test_with_discretized_line(const MBE mbe, const double *p1, const double *p2, double step, int dim, POLY_REGION region, bool equal);
static int region_test_with_discretized_line_count(const MBE mbe, const double *p1, const double *p2, double step, int dim, POLY_REGION region);
static void convert_geom_to_mbe(LWGEOM *geom, MBE *mbe);
static LWGEOM *convert_mbe_to_geom(const SpatialApproximation *ap);
static bool mbe_check_predicate(const SpatialApproximation *ap, const SpatialApproximation *bbox, uint8_t predicate);
static size_t serialize_mbe(const SpatialApproximation *ap, uint8_t *buf);
static SpatialApproximation *recreate_mbe(int n, ApproxPoint *p, uint64_t id);
static void destroy_mbe(SpatialApproximation *ap);

/* ********************
*** Helper functions **
******************** */

double **discretize_line(int *n, const double *p1, const double *p2, double step, int dim)
{
    /*  *n = return matrix row count
        p1/p2 = points (2D double coordinate array)
        step = increment per point
        dim = dimension to increment, x = 0, y = 1
    */
    double **ret;
    double val; // increment value to generate points

    int size = fabs(p1[dim] - p2[dim]) / step, i; //get point count by divinding the absolute distance by the step
    int o_dim = 1;                                //opposite dimension to dim
    if (dim == 1)
    {
        o_dim = 0;
    }
    ret = (double **)lwalloc(size * sizeof(double *)); //matrix of points: n by 2
    if (DB_LT(p1[dim], p2[dim]))                       //incrementing only the smaller point in dim
    {
        val = p1[dim];
    }
    else
    {
        val = p2[dim];
    }
    for (i = 0; i < size; i++) //generating "size" number of points
    {
        val += step;                                    //incrementing
        ret[i] = (double *)lwalloc(2 * sizeof(double)); //allocating array with size 2
        ret[i][dim] = val;                              // setting val to dim dimension
        ret[i][o_dim] = p1[o_dim];                      // o_dim stays the same in all points
    }
    *n = size;  //returning thh size
    return ret; //returning the matrix
}

void free_point_matrix(double **m, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        lwfree(m[i]);
    }
    lwfree(m);
}

bool region_test_with_discretized_line(const MBE mbe, const double *p1, const double *p2, double step, int dim, POLY_REGION region, bool equal)
{ // checks if all discrete points are equal or not to the disired region, when any points if found not to be it return false
    int n, i;
    ApproxPoint p; //ApproxPoint var to use for tests
    double **m = discretize_line(&n, p1, p2, step, dim);
    if (equal) //if all are equal = true, if at least one is not = false
    {
        for (i = 0; i < n; i++)
        {
            p.x = m[i][0]; //building point
            p.y = m[i][1];
            if (mbe_region(mbe, p) != region) //region check
            {
                free_point_matrix(m, n);
                return false; //exits lentghy processing if unecessary
            }
        }
    }
    else //if all are different = true, if at least one is not = false
    {
        for (i = 0; i < n; i++)
        {
            p.x = m[i][0];
            p.y = m[i][1];
            if (mbe_region(mbe, p) == region)
            {
                free_point_matrix(m, n);
                return false;
            }
        }
    }
    free_point_matrix(m, n);
    return true;
}

int region_test_with_discretized_line_count(const MBE mbe, const double *p1, const double *p2, double step, int dim, POLY_REGION region)
{ //counting how many points match with region
    int n, i, count = 0;
    ApproxPoint p;
    double **m = discretize_line(&n, p1, p2, step, dim);
    for (i = 0; i < n; i++)
    {
        p.x = m[i][0];
        p.y = m[i][1];
        if (mbe_region(mbe, p) == region)
        {
            count++;
        }
    }
    free_point_matrix(m, n);
    return count;
}

// MBE create_mbe(int n, ApproxPoint *P)
// {
//     return mbe_create(n, P);
// }

/* ********************
****** Functions ******
******************** */

void convert_geom_to_mbe(LWGEOM *geom, MBE *mbe)
{
    if (geom == NULL)
    {
        *mbe = NULL;
    }
    else
    {
        uint32_t numpoints = lwgeom_count_vertices(geom) - 1, i;
        LWPOINTITERATOR *it = lwpointiterator_create(geom);
        ApproxPoint *parray = lwalloc(numpoints * sizeof(ApproxPoint));
        POINT4D p;
        for (i = 0; i < numpoints; i++)
        {
            if (!lwpointiterator_next(it, &p))
            {
                lwpointiterator_destroy(it);
                lwfree(parray);
                return;
            }
            parray[i].x = p.x;
            parray[i].y = p.y;
        }
        lwpointiterator_destroy(it);
        *mbe = mbe_calculate(numpoints, parray);
        lwfree(parray);
    }
}

LWGEOM *convert_mbe_to_geom(const SpatialApproximation *ap)
{
    int num;
    ApproxPoint *P;
    MBE_APPROX *approx = (void *)ap;
    MBE mbe = approx->mbe;

    if (!mbe)
    {
        return NULL;
    }
    P = mbe_get_points(&num, mbe);

    if (num == 1) //point
    {
        LWPOINT *point = lwpoint_make2d(SRID_UNKNOWN, P[0].x, P[0].y);
        return lwpoint_as_lwgeom(point);
    }
    else if (num == 2) //line
    {
        LWLINE *line;
        POINTARRAY *pa;
        POINT4D pt;
        /* Construct point array */
        pa = ptarray_construct_empty(0, 0, 2);
        /* Assign coordinates to POINT2D array */
        pt.x = P[0].x;
        pt.y = P[0].y;
        ptarray_append_point(pa, &pt, LW_TRUE);
        pt.x = P[1].x;
        pt.y = P[1].y;
        ptarray_append_point(pa, &pt, LW_TRUE);
        /* Construct and serialize linestring */
        line = lwline_construct(SRID_UNKNOWN, NULL, pa);
        return lwline_as_lwgeom(line);
    }
    else //poly
    {
        LWPOLY *poly;
        POINTARRAY *pa = ptarray_construct_empty(0, 0, num + 1);
        POINT4D pt;
        POINTARRAY **ppa = lwalloc(sizeof(POINTARRAY *));
        /* Assign coordinates to point array */
        int i;
        for (i = 0; i < num; i++)
        {
            pt.x = P[i].x;
            pt.y = P[i].y;
            ptarray_append_point(pa, &pt, LW_TRUE);
        }
        //reinserting first point
        pt.x = P[0].x;
        pt.y = P[0].y;
        ptarray_append_point(pa, &pt, LW_TRUE);
        /* Construct polygon */
        ppa[0] = pa;
        poly = lwpoly_construct(SRID_UNKNOWN, NULL, 1, ppa);
        return lwpoly_as_lwgeom(poly);
    }
}

/* do mbe and bbox insersect */
bool intersect(const MBE mbe, const BBox *bbox)
{
    ApproxPoint P1, P2, P3, P4;
    double p1[2], p2[2];
    /*quick test 1: if any of mbe's support points are inside bbox's or on its boundary */
    int n;
    ApproxPoint *P = mbe_get_points(&n, mbe);
    int i = -1;
    do //checking if every point is inside bbox or on it's border
    {
        i++;
        if (DB_LE(bbox->min[0], P[i].x) && DB_GE(bbox->max[0], P[i].x))
        {
            if (DB_LE(bbox->min[1], P[i].y) && DB_GE(bbox->max[1], P[i].y))
            {
                return true;
            }
        }
    } while (i != n);

    /*quick test 2: check if any of bbox's points are inside mbe or on its border */
    P1.x = bbox->min[0];
    P1.y = bbox->min[1];
    P2.x = bbox->min[0];
    P2.y = bbox->max[1];
    P3.x = bbox->max[0];
    P3.y = bbox->max[1];
    P4.x = bbox->max[0];
    P4.y = bbox->min[1];
    if (mbe_region(mbe, P1) != EXTERIOR || mbe_region(mbe, P2) != EXTERIOR || mbe_region(mbe, P3) != EXTERIOR || mbe_region(mbe, P4) != EXTERIOR)
    {
        return true;
    }

    /*slow test: get discrete points for each of bbox's lines and check if anyone is inside mbe or on its border */
    p1[0] = bbox->min[0];
    p1[1] = bbox->max[1];
    p2[0] = bbox->max[0];
    p2[1] = bbox->min[1];

    /*checks for discrete points in bottom, left, right and top lines of bbox, chained conditions with && so it returns in the first true encountered*/
    return region_test_with_discretized_line(mbe, bbox->min, p2, 0.01, 0, EXTERIOR, true) ||
           region_test_with_discretized_line(mbe, bbox->min, p1, 0.01, 1, EXTERIOR, true) ||
           region_test_with_discretized_line(mbe, p2, bbox->max, 0.01, 1, EXTERIOR, true) ||
           region_test_with_discretized_line(mbe, p1, bbox->max, 0.01, 0, EXTERIOR, true);
}

/* mbe is inside OR coveredBy bbox - to check if there is a containment relationship*/
bool inside_or_coveredBy(const MBE mbe, const BBox *bbox)
{
    ApproxPoint P1, P2, P3, P4;
    double p1[2], p2[2];
    /* 1st test - are all the mbe support points inside bbox? */
    int n_mbe;
    ApproxPoint *P = mbe_get_points(&n_mbe, mbe);
    int i = -1;
    do //checking if every point is inside bbox
    {
        i++;
        /* if the point is not inside, then return false */
        if (DB_GE(bbox->min[0], P[i].x) || DB_LE(bbox->max[0], P[i].x))
        {
            if (DB_GE(bbox->min[1], P[i].y) || DB_LE(bbox->max[1], P[i].y))
            {
                return false;
            }
        }
    } while (i != (n_mbe - 1));

    /* 2nd test - are all the bbox points outside mbe or on its border? */
    P1.x = bbox->min[0];
    P1.y = bbox->min[1];
    P2.x = bbox->min[0];
    P2.y = bbox->max[1];
    P3.x = bbox->max[0];
    P3.y = bbox->max[1];
    P4.x = bbox->max[0];
    P4.y = bbox->min[1];
    if (mbe_region(mbe, P1) == INTERIOR && mbe_region(mbe, P2) == INTERIOR && mbe_region(mbe, P3) == INTERIOR && mbe_region(mbe, P4) == INTERIOR)
    {
        return false;
    }

    /*3rd test (slow): are all the discrete bbox points outside mbe or on its border? */
    p1[0] = bbox->min[0];
    p1[1] = bbox->max[1];
    p2[0] = bbox->max[0];
    p2[1] = bbox->min[1];

    /*checks for discrete points in bottom, left, right and top lines of bbox, chained conditions with && so it returns in the first false encountered*/
    return region_test_with_discretized_line(mbe, bbox->min, p2, 0.01, 0, INTERIOR, false) &&
           region_test_with_discretized_line(mbe, bbox->min, p1, 0.01, 1, INTERIOR, false) &&
           region_test_with_discretized_line(mbe, p2, bbox->max, 0.01, 1, INTERIOR, false) &&
           region_test_with_discretized_line(mbe, p1, bbox->max, 0.01, 0, INTERIOR, false);
}

/*is mbe inside bbox?*/
bool inside(const MBE mbe, const BBox *bbox)
{
    ApproxPoint P1, P2, P3, P4;
    double p1[2], p2[2];
    /* 1st test - are all the mbe support points inside bbox? */
    int n_mbe;
    ApproxPoint *P = mbe_get_points(&n_mbe, mbe);
    int i = -1;
    do //checking if every point is inside bbox
    {
        i++;
        /* if the point is not inside, then return false */
        if (DB_GE(bbox->min[0], P[i].x) || DB_LE(bbox->max[0], P[i].x))
        {
            if (DB_GE(bbox->min[1], P[i].y) || DB_LE(bbox->max[1], P[i].y))
            {
                return false;
            }
        }
    } while (i != (n_mbe - 1));

    /* 2nd test - are all the bbox points outside mbe? */
    P1.x = bbox->min[0];
    P1.y = bbox->min[1];
    P2.x = bbox->min[0];
    P2.y = bbox->max[1];
    P3.x = bbox->max[0];
    P3.y = bbox->max[1];
    P4.x = bbox->max[0];
    P4.y = bbox->min[1];
    if (mbe_region(mbe, P1) != EXTERIOR && mbe_region(mbe, P2) != EXTERIOR && mbe_region(mbe, P3) != EXTERIOR && mbe_region(mbe, P4) != EXTERIOR)
    {
        return false;
    }

    /*3rd test (slow): are all the discrete bbox points outside mbe? */
    p1[0] = bbox->min[0];
    p1[1] = bbox->max[1];
    p2[0] = bbox->max[0];
    p2[1] = bbox->min[1];

    /*checks for discrete points in bottom, left, right and top lines of bbox, chained conditions with && so it returns in the first false encountered*/
    return region_test_with_discretized_line(mbe, bbox->min, p2, 0.01, 0, EXTERIOR, true) &&
           region_test_with_discretized_line(mbe, bbox->min, p1, 0.01, 1, EXTERIOR, true) &&
           region_test_with_discretized_line(mbe, p2, bbox->max, 0.01, 1, EXTERIOR, true) &&
           region_test_with_discretized_line(mbe, p1, bbox->max, 0.01, 0, EXTERIOR, true);
}

/*does mbe contains bbox?*/
bool contains(const MBE mbe, const BBox *bbox)
{
    ApproxPoint p1, p2, p3, p4;
    p1.x = bbox->min[0];
    p1.y = bbox->min[1];
    p2.x = bbox->min[0];
    p2.y = bbox->max[1];
    p3.x = bbox->max[0];
    p3.y = bbox->max[1];
    p4.x = bbox->max[0];
    p4.y = bbox->min[1];
    /* are all bbox points inside mbe? */
    if (mbe_is_in_interior(mbe, p1) && mbe_is_in_interior(mbe, p2) && mbe_is_in_interior(mbe, p3) && mbe_is_in_interior(mbe, p4))
    {
        return true;
    }
    return false;
}

/*is mbe covered by bbox?*/
bool coveredBy(const MBE mbe, const BBox *bbox)
{
    int acum = 0;
    /*1st test: any border points are equal? */
    double p1[2], p2[2];
    p1[0] = bbox->min[0];
    p1[1] = bbox->max[1];
    p2[0] = bbox->max[0];
    p2[1] = bbox->min[1];
    acum += region_test_with_discretized_line_count(mbe, bbox->min, p2, 0.01, 0, BOUNDARY); //bottom
    acum += region_test_with_discretized_line_count(mbe, bbox->min, p1, 0.01, 1, BOUNDARY); //left
    acum += region_test_with_discretized_line_count(mbe, p2, bbox->max, 0.01, 1, BOUNDARY); //right
    acum += region_test_with_discretized_line_count(mbe, p1, bbox->max, 0.01, 0, BOUNDARY); //top
    if (acum == 0)
    {
        return false;
    }

    return inside_or_coveredBy(mbe, bbox);
}

/*does mbe cover bbox?*/
bool covers(const MBE mbe, const BBox *bbox)
{
    int borderCount = 0;

    ApproxPoint p1, p2, p3, p4;
    p1.x = bbox->min[0];
    p1.y = bbox->min[1];
    p2.x = bbox->min[0];
    p2.y = bbox->max[1];
    p3.x = bbox->max[0];
    p3.y = bbox->max[1];
    p4.x = bbox->max[0];
    p4.y = bbox->min[1];
    //checking borders
    if (mbe_is_in_boundary(mbe, p1))
    {
        borderCount++;
    }
    if (mbe_is_in_boundary(mbe, p2))
    {
        borderCount++;
    }
    if (mbe_is_in_boundary(mbe, p3))
    {
        borderCount++;
    }
    if (mbe_is_in_boundary(mbe, p4))
    {
        borderCount++;
    }

    if (borderCount == 0)
    {
        return false; //no point is on the border
    }

    /* are all bbox points inside mbe? and do they share any border points? */
    if (mbe_region(mbe, p1) != EXTERIOR && mbe_region(mbe, p2) != EXTERIOR && mbe_region(mbe, p3) != EXTERIOR && mbe_region(mbe, p4) != EXTERIOR)
    {
        return true;
    }
    return false;
}

/* this overlap is according to the 9-IM! */
bool overlap(const MBE mbe, const BBox *bbox)
{
    return !inside_or_coveredBy(mbe, bbox) && !contain_or_covers(mbe, bbox) && intersect(mbe, bbox) && !meet(mbe, bbox);
}

/* this meet is according to the 9-IM! */
bool meet(const MBE mbe, const BBox *bbox)
{
    double p1[2], p2[2];

    /*pre-test: are mbe and bbox not inside one another?*/
    if (inside_or_coveredBy(mbe, bbox) || contain_or_covers(mbe, bbox))
    {
        return false;
    }

    /*1st test (slow): are all the discrete bbox points outside mbe? */
    p1[0] = bbox->min[0];
    p1[1] = bbox->max[1];
    p2[0] = bbox->max[0];
    p2[1] = bbox->min[1];

    /* bottom line */
    if (!region_test_with_discretized_line(mbe, bbox->min, p2, 0.01, 0, INTERIOR, false))
    {
        return false;
    }

    /* left line */
    if (!region_test_with_discretized_line(mbe, bbox->min, p1, 0.01, 1, INTERIOR, false))
    {
        return false;
    }

    /* right line */
    if (!region_test_with_discretized_line(mbe, p2, bbox->max, 0.01, 1, INTERIOR, false))
    {
        return false;
    }

    /* top line */
    if (!region_test_with_discretized_line(mbe, p1, bbox->max, 0.01, 0, INTERIOR, false))
    {
        return false;
    }

    /* if they are not inside or covered by one another and they cannot possibly overlap, if they intersect they meet */
    return intersect(mbe, bbox);
}

/* are mbe and bbox equal? only if both are points */
bool equal(const MBE mbe, const BBox *bbox)
{
    /*check if mbe number of points is 1 */
    int num;
    ApproxPoint *P = mbe_get_points(&num, mbe);
    if (num == 1)
    {
        /*check if bbox is a point */
        if (DB_IS_EQUAL(bbox->min[0], bbox->max[0]) && DB_IS_EQUAL(bbox->min[1], bbox->max[1]))
        {
            /* check if points match */
            if (DB_IS_EQUAL(bbox->min[0], P[0].x) && DB_IS_EQUAL(bbox->min[1], P[0].y))
            {
                return true;
            }
        }
    }
    return false;
}

bool contain_or_covers(const MBE mbe, const BBox *bbox)
{
    ApproxPoint p1, p2, p3, p4;
    p1.x = bbox->min[0];
    p1.y = bbox->min[1];
    p2.x = bbox->min[0];
    p2.y = bbox->max[1];
    p3.x = bbox->max[0];
    p3.y = bbox->max[1];
    p4.x = bbox->max[0];
    p4.y = bbox->min[1];
    /* are all bbox points inside mbe? or do they share any border ponits? */
    if (mbe_region(mbe, p1) != EXTERIOR && mbe_region(mbe, p2) != EXTERIOR && mbe_region(mbe, p3) != EXTERIOR && mbe_region(mbe, p4) != EXTERIOR)
    {
        return true;
    }
    return false;
}

bool mbe_check_predicate(const SpatialApproximation *ap, const SpatialApproximation *bbox, uint8_t predicate)
{
    if (bbox->type == BBOX_TYPE)
    {
        MBE_APPROX *mbe = (void *)ap;
        BBOX_APPROX *bbox_ap = (void *)bbox;
        BBox *bbox1 = &(bbox_ap->bbox);

        switch (predicate)
        {
        case INTERSECTS:
            return intersect(mbe->mbe, bbox1);
        case DISJOINT:
            return !intersect(mbe->mbe, bbox1);
        case OVERLAP:
            return overlap(mbe->mbe, bbox1);
        case MEET:
            return meet(mbe->mbe, bbox1);
        case INSIDE:
            return inside(mbe->mbe, bbox1);
        case CONTAINS:
            return contains(mbe->mbe, bbox1);
        case COVEREDBY:
            return coveredBy(mbe->mbe, bbox1);
        case COVERS:
            return covers(mbe->mbe, bbox1);
        case EQUAL:
            return equal(mbe->mbe, bbox1);
        case INSIDE_OR_COVEREDBY:
            return inside_or_coveredBy(mbe->mbe, bbox1);
        case CONTAINS_OR_COVERS:
            return contain_or_covers(mbe->mbe, bbox1);
        default:
            return false;
        }
        return false;
    }
    else
    {
        _DEBUG(WARNING, "Wrong type for bbox approximation\n");
        return false;
    }
}

size_t serialize_mbe(const SpatialApproximation *ap, uint8_t *buf)
{
    ApproxPoint nil = {0, 0};
    ApproxPoint *p;
    uint8_t *loc, u_n;
    int i;
    int n;
    MBE_APPROX *mbe_ap = (void *)ap;
    loc = buf;

    memcpy(loc, &(ap->id), sizeof(uint64_t));
    loc += sizeof(uint64_t);

    p = mbe_get_points(&n, mbe_ap->mbe);

    u_n = n;
    memcpy(loc, &u_n, sizeof(uint8_t));
    loc += sizeof(uint8_t);

    for (i = 0; i < n; i++)
    {
        memcpy(loc, &(p[i]), sizeof(ApproxPoint));
        loc += sizeof(ApproxPoint);
    }

    for (i = 0; i < (5 - n); i++)
    {
        memcpy(loc, &nil, sizeof(ApproxPoint));
        loc += sizeof(ApproxPoint);
    }

    return (size_t)(loc - buf);
}

size_t mbe_size()
{
    return (sizeof(ApproxPoint) * 5) + sizeof(uint64_t) + sizeof(uint8_t);
}

bool write_mbe_to_file(const FileSpecification *fs, const SpatialApproximation **ap, int n)
{
    uint8_t *buf, *loc;
    int i;
    size_t total_size = (n * mbe_size()) + sizeof(int), finalsize;

    if (total_size > fs->page_size)
    {
        _DEBUG(ERROR, "Total size is bigger than page size\n");
        return false;
    }

    if (fs->io_access == DIRECT_ACCESS)
    {
        if (posix_memalign((void **)&buf, fs->page_size, fs->page_size))
        {
            _DEBUG(ERROR, "Allocation failed at MBE write_to_file\n");
            return false;
        }
    }
    else
    {
        buf = lwalloc(fs->page_size);
    }

    loc = buf;

    memcpy(loc, &n, sizeof(int));
    loc += sizeof(int);

    for (i = 0; i < n; i++)
    {
        loc += serialize_mbe(ap[i], loc);
    }

    finalsize = (size_t)(loc - buf);

    if (total_size != finalsize)
    {
        _DEBUGF(WARNING, "Return size (%zu) not equal to expected size (%zu)\n", finalsize, total_size);
        return false;
    }
    append_page(fs, buf);

    if (fs->io_access == DIRECT_ACCESS)
    {
        free(buf);
    }
    else
    {
        lwfree(buf);
    }

    return true;
}

SpatialApproximation *recreate_mbe(int n, ApproxPoint *p, uint64_t id)
{
    static const SpatialApproximationInterface vtable = {convert_mbe_to_geom, mbe_check_predicate, destroy_mbe};
    static SpatialApproximation base = {&vtable};
    MBE_APPROX *approx;

    MBE mbe;
    mbe = mbe_create(n, p);

    base.type = MBE_TYPE;
    base.id = id;
    approx = lwalloc(sizeof(MBE_APPROX));
    memcpy(&approx->base, &base, sizeof(base));
    approx->mbe = mbe;

    return &approx->base;
}

SpatialApproximation *read_serialized_mbe(uint8_t *buf, size_t *size)
{
    SpatialApproximation *ap;
    uint8_t *start_ptr = buf;
    uint8_t n;
    uint64_t id;
    ApproxPoint p[5];
    int i;

    memcpy(&id, buf, sizeof(uint64_t));
    buf += sizeof(uint64_t);

    memcpy(&n, buf, sizeof(uint8_t));
    buf += sizeof(uint8_t);

    for (i = 0; i < n; i++)
    {
        memcpy(&(p[i]), buf, sizeof(ApproxPoint));
        buf += sizeof(ApproxPoint);
    }

    buf += sizeof(ApproxPoint) * (5 - n);

    ap = recreate_mbe(n, p, id);

    *size = (size_t)(buf - start_ptr);
    return ap;
}

SpatialApproximation **read_mbe_from_file(const FileSpecification *fs, int page, int *n)
{
    SpatialApproximation **ap;
    uint8_t *buf, *loc;
    int i;
    size_t size;

    if (fs->io_access == DIRECT_ACCESS)
    {
        if (posix_memalign((void **)&buf, fs->page_size, fs->page_size))
        {
            _DEBUG(ERROR, "Allocation failed at MBE read_from_file\n");
            return NULL;
        }
    }
    else
    {
        buf = lwalloc(fs->page_size);
    }

    disk_read_one_page(fs, page, buf);

    loc = buf;

    memcpy(n, loc, sizeof(int));
    loc += sizeof(int);

    if (*n <= 0)
    {
        _DEBUGF(WARNING, "Number of approximations in page is invalid: %d\n", *n);
        return NULL;
    }

    ap = lwalloc((*n) * sizeof(SpatialApproximation *));

    for (i = 0; i < (*n); i++)
    {
        ap[i] = read_serialized_mbe(loc, &size);
        loc += size;
    }

    if (fs->io_access == DIRECT_ACCESS)
    {
        free(buf);
    }
    else
    {
        lwfree(buf);
    }

    return ap;
}

SpatialApproximation *create_mbe(LWGEOM *geom, uint64_t id)
{
    static const SpatialApproximationInterface vtable = {convert_mbe_to_geom, mbe_check_predicate, destroy_mbe};
    static SpatialApproximation base = {&vtable};
    MBE_APPROX *approx;

    MBE mbe = NULL;
    convert_geom_to_mbe(geom, &mbe);

    base.type = MBE_TYPE;
    base.id = id;
    approx = lwalloc(sizeof(MBE_APPROX));
    memcpy(&approx->base, &base, sizeof(base));
    approx->mbe = mbe;

    return &approx->base;
}

void destroy_mbe(SpatialApproximation *ap)
{
    MBE_APPROX *mbe_ap = (void *)ap;
    mbe_free(mbe_ap->mbe);
    lwfree(mbe_ap);
}