#include <string.h>
#include "../main/lwgeom_handler.h"
#include "rmbp_handler.h"
#include "../main/math_util.h"
#include "unistd.h"
#include "../main/log_messages.h"

static bool intersect(const RMBP rmbp, const BBox *bbox);
static bool overlap(const RMBP rmbp, const BBox *bbox);
static bool meet(const RMBP rmbp, const BBox *bbox);
static bool coveredBy(const RMBP rmbp, const BBox *bbox);
static bool inside(const RMBP rmbp, const BBox *bbox);
static bool covers(const RMBP rmbp, const BBox *bbox);
static bool contains(const RMBP rmbp, const BBox *bbox);
static bool inside_or_coveredBy(const RMBP rmbp, const BBox *bbox);
static bool contain_or_covers(const RMBP rmbp, const BBox *bbox);
static bool equal(const RMBP rmbp, const BBox *bbox);
static void convert_geom_to_rmbp(LWGEOM *geom, RMBP *rmbp);
static LWGEOM *convert_rmbp_to_geom(const SpatialApproximation *ap);
static bool rmbp_check_predicate(const SpatialApproximation *ap, const SpatialApproximation *bbox, uint8_t predicate);
static size_t serialize_rmbp(const SpatialApproximation *ap, uint8_t *buf);
static SpatialApproximation *recreate_rmbp(ApproxPoint p1, ApproxPoint p2, ApproxPoint p3, ApproxPoint p4, uint64_t id);
static void destroy_rmbp(SpatialApproximation *ap);

void convert_geom_to_rmbp(LWGEOM *geom, RMBP *rmbp)
{
    if (geom == NULL)
    {
        *rmbp = NULL;
    }
    else
    {
        uint32_t numpoints, i;
        LWPOINTITERATOR *it;
        ApproxPoint *parray;
        POINT4D p;

        LWGEOM *chgeom = convexhull(geom);
        lwgeom_force_counterclockwise(chgeom);
        numpoints = lwgeom_count_vertices(chgeom) - 1;
        it = lwpointiterator_create(chgeom);
        parray = lwalloc(numpoints * sizeof(ApproxPoint));

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
        *rmbp = rmbp_calculate(numpoints, parray);
        lwfree(parray);
        lwgeom_free(chgeom);
    }
}

LWGEOM *convert_rmbp_to_geom(const SpatialApproximation *ap)
{
    RMBP_APPROX *approx = (void *)ap;
    RMBP rmbp = approx->rmbp;
    ApproxPoint *P;
    double pos;

    if (!rmbp)
    {
        return NULL;
    }

    P = rmbp_get_points(rmbp);
    pos = P[0].x;

    if ((P[0].y == pos) && (P[1].x == pos) && (P[1].y == pos) && (P[2].x == pos) && (P[2].y == pos) && (P[3].x == pos) && (P[3].y == pos)) //point
    {
        LWPOINT *point = lwpoint_make2d(SRID_UNKNOWN, pos, pos);
        return lwpoint_as_lwgeom(point);
    }
    else if ((P[0].x == P[2].x && P[0].y == P[2].y) && (P[1].x == P[3].x && P[1].y == P[3].y)) //line
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
        POINTARRAY *pa = ptarray_construct_empty(0, 0, 5);
        POINT4D pt;
        POINTARRAY **ppa = lwalloc(sizeof(POINTARRAY *));
        /* Assign coordinates to point array */
        pt.x = P[0].x;
        pt.y = P[0].y;
        ptarray_append_point(pa, &pt, LW_TRUE);
        pt.x = P[1].x;
        pt.y = P[1].y;
        ptarray_append_point(pa, &pt, LW_TRUE);
        pt.x = P[2].x;
        pt.y = P[2].y;
        ptarray_append_point(pa, &pt, LW_TRUE);
        pt.x = P[3].x;
        pt.y = P[3].y;
        ptarray_append_point(pa, &pt, LW_TRUE);
        pt.x = P[0].x;
        pt.y = P[0].y;
        ptarray_append_point(pa, &pt, LW_TRUE);
        /* Construct polygon */
        ppa[0] = pa;
        poly = lwpoly_construct(SRID_UNKNOWN, NULL, 1, ppa);
        return lwpoly_as_lwgeom(poly);
    }
}

bool intersect(const RMBP rmbp, const BBox *bbox)
{
    ApproxPoint *P;
    ApproxPoint p1, p2, p3, p4;
    int i = -1;
    /* test 0. check if bbox "center" is inside rmbp */
    ApproxPoint p;
    p.x = bbox->min[0] + (bbox->max[0] - bbox->min[0]) / 2.0;
    p.y = bbox->min[1] + (bbox->max[1] - bbox->min[1]) / 2.0;
    if (rmbp_region(rmbp, p) != EXTERIOR)
    {
        return true;
    }

    /* 1. check if rmbp has any point within bbox boundary */
    P = rmbp_get_points(rmbp);
    do
    {
        i++;
        /* if they share any point, then return true */
        if (DB_LE(bbox->min[0], P[i].x) && DB_GE(bbox->max[0], P[i].x))
        {
            if (DB_LE(bbox->min[1], P[i].y) && DB_GE(bbox->max[1], P[i].y))
            {
                return true;
            }
        }
    } while (i != 3);

    /* 2. check if bbox has any point within rmbp boundary */
    p1.x = bbox->min[0];
    p1.y = bbox->min[1];
    p2.x = bbox->min[0];
    p2.y = bbox->max[1];
    p3.x = bbox->max[0];
    p3.y = bbox->max[1];
    p4.x = bbox->max[0];
    p4.y = bbox->min[1];
    if ((rmbp_region(rmbp, p1) != EXTERIOR) || (rmbp_region(rmbp, p2) != EXTERIOR) || (rmbp_region(rmbp, p3) != EXTERIOR) || (rmbp_region(rmbp, p4) != EXTERIOR))
    {
        return true;
    }
    return false;
}

/* rmbp is inside OR coveredBy bbox - to check if there is a containment relationship*/
bool inside_or_coveredBy(const RMBP rmbp, const BBox *bbox)
{
    ApproxPoint *P = rmbp_get_points(rmbp);
    int i = -1;
    do //checking if every point is inside bbox or on it's border
    {
        i++;
        /* if the point is not inside nor on the border, then return false */
        if (DB_GT(bbox->min[0], P[i].x) || DB_LT(bbox->max[0], P[i].x))
        {
            return false;
        }
        if (DB_GT(bbox->min[1], P[i].y) || DB_LT(bbox->max[1], P[i].y))
        {
            return false;
        }
    } while (i != 3);
    return !equal(rmbp, bbox);
}

/*is rmbp inside bbox?*/
bool inside(const RMBP rmbp, const BBox *bbox)
{
    ApproxPoint *P = rmbp_get_points(rmbp);
    int i = -1;
    do //checking if every point is inside bbox
    {
        i++;
        /* if the point is not inside, then return false */
        if (DB_GE(bbox->min[0], P[i].x) || DB_LE(bbox->max[0], P[i].x))
        {
            return false;
        }
        if (DB_GE(bbox->min[1], P[i].y) || DB_LE(bbox->max[1], P[i].y))
        {
            return false;
        }
    } while (i != 3);
    return true;
}

/*does rmbp contains bbox?*/
bool contains(const RMBP rmbp, const BBox *bbox)
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
    /* are all bbox points inside rmbp? */
    if (rmbp_is_in_interior(rmbp, p1) && rmbp_is_in_interior(rmbp, p2) && rmbp_is_in_interior(rmbp, p3) && rmbp_is_in_interior(rmbp, p4))
    {
        return true;
    }
    return false;
}

/*is rmbp covered by bbox?*/
bool coveredBy(const RMBP rmbp, const BBox *bbox)
{
    ApproxPoint *P = rmbp_get_points(rmbp);
    int i = -1;
    int borderCount = 0;
    do //checking if every point is inside bbox or on it's border, having at least one on the border
    {
        i++;
        /* if the point is not inside nor on the border, then return false */
        if (DB_GT(bbox->min[0], P[i].x) || DB_LT(bbox->max[0], P[i].x))
        {
            return false;
        }
        if (DB_GT(bbox->min[1], P[i].y) || DB_LT(bbox->max[1], P[i].y))
        {
            return false;
        }
        if (DB_IS_EQUAL(bbox->min[0], P[i].x) || DB_IS_EQUAL(bbox->max[0], P[i].x) || DB_IS_EQUAL(bbox->min[1], P[i].y) || DB_IS_EQUAL(bbox->max[1], P[i].y))
        {
            borderCount++;
        }
    } while (i != 3);

    if (borderCount == 0 || (borderCount == 4 && equal(rmbp, bbox)))
    {
        return false;
    }
    else
    {
        return true;
    }
}

/*does rmbp cover bbox?*/
bool covers(const RMBP rmbp, const BBox *bbox)
{
    ApproxPoint p1, p2, p3, p4;
    int borderCount = 0;

    p1.x = bbox->min[0];
    p1.y = bbox->min[1];
    p2.x = bbox->min[0];
    p2.y = bbox->max[1];
    p3.x = bbox->max[0];
    p3.y = bbox->max[1];
    p4.x = bbox->max[0];
    p4.y = bbox->min[1];

    //checking borders
    if (rmbp_is_in_boundary(rmbp, p1))
    {
        borderCount++;
    }
    if (rmbp_is_in_boundary(rmbp, p2))
    {
        borderCount++;
    }
    if (rmbp_is_in_boundary(rmbp, p3))
    {
        borderCount++;
    }
    if (rmbp_is_in_boundary(rmbp, p4))
    {
        borderCount++;
    }

    if (borderCount == 0)
    {
        return false; //no point is on the border
    }
    if (borderCount == 4 && equal(rmbp, bbox))
    {
        return false; //the rectangles are the same
    }

    /* are all bbox points inside rmbp? and do they share any border points? */
    if ((rmbp_region(rmbp, p1) != EXTERIOR) && (rmbp_region(rmbp, p2) != EXTERIOR) && (rmbp_region(rmbp, p3) != EXTERIOR) && (rmbp_region(rmbp, p4) != EXTERIOR))
    {
        return true;
    }
    return false;
}

/* this overlap is according to the 9-IM! */
bool overlap(const RMBP rmbp, const BBox *bbox)
{
    return !inside_or_coveredBy(rmbp, bbox) && !contain_or_covers(rmbp, bbox) && intersect(rmbp, bbox) && !meet(rmbp, bbox) && !equal(rmbp, bbox);
}

/* this meet is according to the 9-IM! */
bool meet(const RMBP rmbp, const BBox *bbox)
{
    ApproxPoint *P = rmbp_get_points(rmbp);
    int c1 = 0, c2 = 0, c3 = 0;
    ApproxPoint p;
    int i = -1;
    ApproxPoint p1, p2, p3, p4;

    /*pre-test: are rmbp and bbox not inside one another?*/
    if (inside_or_coveredBy(rmbp, bbox) || contain_or_covers(rmbp, bbox))
    {
        return false;
    }

    /* test 0. check if bbox "center" is inside rmbp */
    p.x = bbox->min[0] + (bbox->max[0] - bbox->min[0]) / 2;
    p.y = bbox->min[1] + (bbox->max[1] - bbox->min[1]) / 2;
    if (rmbp_region(rmbp, p) != EXTERIOR)
    {
        return false;
    }

    /*1st test - validate if rmbp touches bbox on the Y axis*/
    do
    {
        i++;
        /* if the point is equal to xmin or xmax and is in the range of y */
        if ((DB_IS_EQUAL(bbox->min[0], P[i].x) || DB_IS_EQUAL(bbox->max[0], P[i].x)) && (DB_LE(bbox->min[1], P[i].y) && DB_GE(bbox->max[1], P[i].y)))
        {
            c1++;
        }
    } while (i != 3);
    if (c1 > 2) //rmbp touches bbox in too many points (overlap/coveredBy/intersects/inside)
    {
        return false;
    }
    else if (c1 > 0) //rmbp touches bbox in one or two points in the x axis
    {
        return true;
    }

    /*2nd test - validate if rmbp touches bbox on the X axis*/
    i = -1;
    do
    {
        i++;
        /* if the point is equal to ymin or ymax and is in the range of x*/
        if ((DB_IS_EQUAL(bbox->min[1], P[i].y) || DB_IS_EQUAL(bbox->max[1], P[i].y)) && (DB_LE(bbox->min[0], P[i].x) && DB_GE(bbox->max[0], P[i].x)))
        {
            c2++;
        }
    } while (i != 3);
    if (c2 > 2) //rmbp touches bbox in too many points (overlap/coveredBy/intersects/inside)
    {
        return false;
    }
    else if (c2 > 0) //rmbp touches bbox in one or two points in the y axis
    {
        return true;
    }

    /*3rd test - validate if bbox touches rmbp borders*/
    p1.x = bbox->min[0];
    p1.y = bbox->min[1];
    p2.x = bbox->min[0];
    p2.y = bbox->max[1];
    p3.x = bbox->max[0];
    p3.y = bbox->max[1];
    p4.x = bbox->max[0];
    p4.y = bbox->min[1];

    /* check if bbox point is in rmbp border*/
    if (rmbp_is_in_boundary(rmbp, p1))
    {
        c3++;
    }
    if (rmbp_is_in_boundary(rmbp, p2))
    {
        c3++;
    }
    if (rmbp_is_in_boundary(rmbp, p3))
    {
        c3++;
    }
    if (rmbp_is_in_boundary(rmbp, p4))
    {
        c3++;
    }

    if (c3 == 0) //rmbp is not touched by bbox in it's border
    {
        return false;
    }
    else if (c3 > 2) //bbox touches rmbp in too many points (overlap/coveredBy/intersects/inside)
    {
        return false;
    }
    else //bbox touches rmbp border in one or two points
    {
        return true;
    }
}

bool equal(const RMBP rmbp, const BBox *bbox)
{
    ApproxPoint *P = rmbp_get_points(rmbp);
    int i, u1 = -1, u2 = -1, u3 = -1; //-1 equals not select
    //check bbox point 1
    ApproxPoint p1, p2, p3, p4;
    p1.x = bbox->min[0];
    p1.y = bbox->min[1];
    for (i = 0; i < 4; i++)
    {
        if (DB_IS_EQUAL(P[i].x, p1.x) && DB_IS_EQUAL(P[i].y, p1.y))
        {
            u1 = i;
            break;
        }
    }
    if (u1 == -1) //if u1 is not selected no point is the same as 1, so not equal
    {
        return false;
    }
    //check bbox point 2
    p2.x = bbox->min[0];
    p2.y = bbox->max[1];
    for (int i = 0; i < 4; i++)
    {
        if (i != u1)
        {
            if (DB_IS_EQUAL(P[i].x, p2.x) && DB_IS_EQUAL(P[i].y, p2.y))
            {
                u2 = i;
                break;
            }
        }
    }
    if (u2 == -1) //if u2 is not selected no point is the same as 2, so not equal
    {
        return false;
    }
    //check bbox point 3
    p3.x = bbox->max[0];
    p3.y = bbox->max[1];
    for (int i = 0; i < 4; i++)
    {
        if (i != u1 && i != u2)
        {
            if (DB_IS_EQUAL(P[i].x, p3.x) && DB_IS_EQUAL(P[i].y, p3.y))
            {
                u3 = i;
                break;
            }
        }
    }
    if (u3 == -1) //if u3 is not selected no point is the same as 3, so not equal
    {
        return false;
    }
    //check bbox point 4
    p4.x = bbox->max[0];
    p4.y = bbox->min[1];
    for (int i = 0; i < 4; i++)
    {
        if (i != u1 && i != u2 && i != u3)
        {
            if (DB_IS_EQUAL(P[i].x, p4.x) && DB_IS_EQUAL(P[i].y, p4.y))
            {
                return true; //if the last rmbp point matches the last bbox point
            }
            else
            {
                return false; //the last points don't match
            }
        }
    }
    return false;
}

bool contain_or_covers(const RMBP rmbp, const BBox *bbox)
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
    /* are all bbox points inside rmbp? or do they share any border points? */
    if ((rmbp_region(rmbp, p1) != EXTERIOR) && (rmbp_region(rmbp, p2) != EXTERIOR) && (rmbp_region(rmbp, p3) != EXTERIOR) && (rmbp_region(rmbp, p4) != EXTERIOR))
    {
        return !equal(rmbp, bbox);
    }
    return false;
}

bool rmbp_check_predicate(const SpatialApproximation *ap, const SpatialApproximation *bbox, uint8_t predicate)
{
    if (bbox->type == BBOX_TYPE)
    {
        RMBP_APPROX *rmbp = (void *)ap;
        BBOX_APPROX *bbox_ap = (void *)bbox;
        BBox *bbox1 = &(bbox_ap->bbox);
        switch (predicate)
        {
        case INTERSECTS:
            return intersect(rmbp->rmbp, bbox1);
        case DISJOINT:
            return !intersect(rmbp->rmbp, bbox1);
        case OVERLAP:
            return overlap(rmbp->rmbp, bbox1);
        case MEET:
            return meet(rmbp->rmbp, bbox1);
        case INSIDE:
            return inside(rmbp->rmbp, bbox1);
        case CONTAINS:
            return contains(rmbp->rmbp, bbox1);
        case COVEREDBY:
            return coveredBy(rmbp->rmbp, bbox1);
        case COVERS:
            return covers(rmbp->rmbp, bbox1);
        case EQUAL:
            return equal(rmbp->rmbp, bbox1);
        case INSIDE_OR_COVEREDBY:
            return inside_or_coveredBy(rmbp->rmbp, bbox1);
        case CONTAINS_OR_COVERS:
            return contain_or_covers(rmbp->rmbp, bbox1);
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

size_t serialize_rmbp(const SpatialApproximation *ap, uint8_t *buf)
{
    uint8_t *loc;
    int i;
    ApproxPoint *p;
    RMBP_APPROX *rmbp_ap = (void *)ap;
    loc = buf;

    memcpy(loc, &(ap->id), sizeof(uint64_t));
    loc += sizeof(uint64_t);

    p = rmbp_get_points(rmbp_ap->rmbp);

    for (i = 0; i < 4; i++)
    {
        memcpy(loc, &(p[i]), sizeof(ApproxPoint));
        loc += sizeof(ApproxPoint);
    }
    return (size_t)(loc - buf);
}

size_t rmbp_size()
{
    return (sizeof(ApproxPoint) * 4) + sizeof(uint64_t);
}

bool write_rmbp_to_file(const FileSpecification *fs, const SpatialApproximation **ap, int n)
{
    uint8_t *buf, *loc;
    int i;
    size_t total_size = (n * rmbp_size()) + sizeof(int), finalsize;

    if (total_size > fs->page_size)
    {
        _DEBUG(ERROR, "Total size is bigger than page size\n");
        return false;
    }

    if (fs->io_access == DIRECT_ACCESS)
    {
        if (posix_memalign((void **)&buf, fs->page_size, fs->page_size))
        {
            _DEBUG(ERROR, "Allocation failed at RMBP write_to_file\n");
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
        loc += serialize_rmbp(ap[i], loc);
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

SpatialApproximation *recreate_rmbp(ApproxPoint p1, ApproxPoint p2, ApproxPoint p3, ApproxPoint p4, uint64_t id)
{
    static const SpatialApproximationInterface vtable = {convert_rmbp_to_geom, rmbp_check_predicate, destroy_rmbp};
    static SpatialApproximation base = {&vtable};
    RMBP_APPROX *approx;

    RMBP rmbp;
    rmbp = rmbp_create(p1, p2, p3, p4);

    base.type = RMBP_TYPE;
    base.id = id;
    approx = lwalloc(sizeof(RMBP_APPROX));
    memcpy(&approx->base, &base, sizeof(base));
    approx->rmbp = rmbp;

    return &approx->base;
}

SpatialApproximation *read_serialized_rmbp(uint8_t *buf, size_t *size)
{
    uint8_t *start_ptr = buf;
    ApproxPoint p[4];
    uint64_t id;
    int i;
    SpatialApproximation *ap;

    memcpy(&id, buf, sizeof(uint64_t));
    buf += sizeof(uint64_t);

    for (i = 0; i < 4; i++)
    {
        memcpy(&(p[i]), buf, sizeof(ApproxPoint));
        buf += sizeof(ApproxPoint);
    }

    ap = recreate_rmbp(p[0], p[1], p[2], p[3], id);

    *size = (size_t)(buf - start_ptr);
    return ap;
}

SpatialApproximation **read_rmbp_from_file(const FileSpecification *fs, int page, int *n)
{
    uint8_t *buf, *loc;
    SpatialApproximation **ap;
    int i;
    size_t size;

    if (fs->io_access == DIRECT_ACCESS)
    {
        if (posix_memalign((void **)&buf, fs->page_size, fs->page_size))
        {
            _DEBUG(ERROR, "Allocation failed at RMBP read_from_file\n");
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
        ap[i] = read_serialized_rmbp(loc, &size);
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

SpatialApproximation *create_rmbp(LWGEOM *geom, uint64_t id)
{
    static const SpatialApproximationInterface vtable = {convert_rmbp_to_geom, rmbp_check_predicate, destroy_rmbp};
    static SpatialApproximation base = {&vtable};
    RMBP_APPROX *approx;

    RMBP rmbp = NULL;
    convert_geom_to_rmbp(geom, &rmbp);

    base.type = RMBP_TYPE;
    base.id = id;
    approx = lwalloc(sizeof(RMBP_APPROX));
    memcpy(&approx->base, &base, sizeof(base));
    approx->rmbp = rmbp;

    return &approx->base;
}

void destroy_rmbp(SpatialApproximation *ap)
{
    RMBP_APPROX *rmbp_ap = (void *)ap;
    rmbp_free(rmbp_ap->rmbp);
    lwfree(rmbp_ap);
}