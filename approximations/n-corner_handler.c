#include <string.h>
#include "lwgeom_log.h"
#include "lwgeom_geos.h"
#include "../main/lwgeom_handler.h"
#include "n-corner_handler.h"
#include "../main/math_util.h"
#include "../main/log_messages.h"

static bool intersect(const N_CORNER *n_corner, const BBox *bbox);
static bool overlap(const N_CORNER *n_corner, const BBox *bbox);
static bool meet(const N_CORNER *n_corner, const BBox *bbox);
static bool coveredBy(const N_CORNER *n_corner, const BBox *bbox);
static bool inside(const N_CORNER *n_corner, const BBox *bbox);
static bool covers(const N_CORNER *n_corner, const BBox *bbox);
static bool contains(const N_CORNER *n_corner, const BBox *bbox);
static bool inside_or_coveredBy(const N_CORNER *n_corner, const BBox *bbox);
static bool contain_or_covers(const N_CORNER *n_corner, const BBox *bbox);
static bool equal(const N_CORNER *n_corner, const BBox *bbox);
static void convert_geom_to_5_corner(LWGEOM *geom, N_CORNER *n_corner);
static void convert_geom_to_4_corner(LWGEOM *geom, N_CORNER *n_corner);
static LWGEOM *convert_n_corner_to_geom(const SpatialApproximation *ap);
static LWGEOM *convert_n_corner_to_geom_internal(const N_CORNER *n_corner);
static bool n_corner_check_predicate(const SpatialApproximation *ap, const SpatialApproximation *bbox, uint8_t predicate);
static size_t serialize_n_corner(const SpatialApproximation *ap, uint8_t *buf);
static size_t n_corner_size(uint8_t type);
static SpatialApproximation *recreate_n_corner(int n, ApproxPoint *p, uint8_t type, uint64_t id);
static void destroy_n_corner(SpatialApproximation *ap);

void convert_geom_to_5_corner(LWGEOM *geom, N_CORNER *n_corner)
{
    if (geom == NULL)
    {
        n_corner = NULL;
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
        *n_corner = five_corner_calculate(numpoints, parray);
        lwfree(parray);
        lwgeom_free(chgeom);
    }
}

void convert_geom_to_4_corner(LWGEOM *geom, N_CORNER *n_corner)
{
    if (geom == NULL)
    {
        n_corner = NULL;
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
        *n_corner = four_corner_calculate(numpoints, parray);
        lwfree(parray);
        lwgeom_free(chgeom);
    }
}

LWGEOM *convert_n_corner_to_geom(const SpatialApproximation *ap)
{
    ApproxPoint *P;
    N_CORNER_APPROX *approx = (void *)ap;
    N_CORNER *n_corner = &(approx->n_corner);

    if (!n_corner)
    {
        return NULL;
    }

    P = n_corner->p;

    if (n_corner->num == 1) //point
    {
        LWPOINT *point = lwpoint_make2d(SRID_UNKNOWN, P[0].x, P[0].y);
        return lwpoint_as_lwgeom(point);
    }
    else if (n_corner->num == 2) //line
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
        POINTARRAY *pa = ptarray_construct_empty(0, 0, (n_corner->num) + 1);
        POINT4D pt;
        POINTARRAY **ppa = lwalloc(sizeof(POINTARRAY *));
        /* Assign coordinates to point array */
        int i;
        for (i = 0; i < (n_corner->num); i++)
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

LWGEOM *convert_n_corner_to_geom_internal(const N_CORNER *n_corner)
{
    ApproxPoint *P;
    if (!n_corner)
    {
        return NULL;
    }

    P = n_corner->p;

    if (n_corner->num == 1) //point
    {
        LWPOINT *point = lwpoint_make2d(SRID_UNKNOWN, P[0].x, P[0].y);
        return lwpoint_as_lwgeom(point);
    }
    else if (n_corner->num == 2) //line
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
        POINTARRAY *pa = ptarray_construct_empty(0, 0, (n_corner->num) + 1);
        POINT4D pt;
        POINTARRAY **ppa = lwalloc(sizeof(POINTARRAY *));
        /* Assign coordinates to point array */
        int i;
        for (i = 0; i < (n_corner->num); i++)
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

bool intersect(const N_CORNER *n_corner, const BBox *bbox)
{
    int result;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    LWGEOM *geom1 = convert_n_corner_to_geom_internal(n_corner);
    LWGEOM *geom2 = convert_bbox_to_geom_internal(bbox);

    initGEOS(lwnotice, lwgeom_geos_error);
    g1 = (GEOSGeometry *)LWGEOM2GEOS(geom1, LW_FALSE);
    g2 = (GEOSGeometry *)LWGEOM2GEOS(geom2, LW_FALSE);

    result = GEOSIntersects(g1, g2);
    GEOSGeom_destroy(g1);
    GEOSGeom_destroy(g2);
    lwgeom_free(geom1);
    lwgeom_free(geom2);
    return result;
}

/* n_corner is inside OR coveredBy bbox - to check if there is a containment relationship*/
bool inside_or_coveredBy(const N_CORNER *n_corner, const BBox *bbox)
{
    int result;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    LWGEOM *geom1 = convert_n_corner_to_geom_internal(n_corner);
    LWGEOM *geom2 = convert_bbox_to_geom_internal(bbox);

    initGEOS(lwnotice, lwgeom_geos_error);
    g1 = (GEOSGeometry *)LWGEOM2GEOS(geom1, LW_FALSE);
    g2 = (GEOSGeometry *)LWGEOM2GEOS(geom2, LW_FALSE);

    result = GEOSWithin(g1, g2);
    GEOSGeom_destroy(g1);
    GEOSGeom_destroy(g2);
    lwgeom_free(geom1);
    lwgeom_free(geom2);
    return result && !equal(n_corner, bbox);
}

/*is n_corner inside bbox?*/
bool inside(const N_CORNER *n_corner, const BBox *bbox)
{
    int result;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    LWGEOM *geom1 = convert_n_corner_to_geom_internal(n_corner);
    LWGEOM *geom2 = convert_bbox_to_geom_internal(bbox);

    initGEOS(lwnotice, lwgeom_geos_error);
    g1 = (GEOSGeometry *)LWGEOM2GEOS(geom1, LW_FALSE);
    g2 = (GEOSGeometry *)LWGEOM2GEOS(geom2, LW_FALSE);

    result = GEOSRelatePattern(g1, g2, "T*F*F*T**");
    GEOSGeom_destroy(g1);
    GEOSGeom_destroy(g2);
    lwgeom_free(geom1);
    lwgeom_free(geom2);
    return result;
}

/*does n_corner contains bbox?*/
bool contains(const N_CORNER *n_corner, const BBox *bbox)
{
    int result;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    LWGEOM *geom1 = convert_n_corner_to_geom_internal(n_corner);
    LWGEOM *geom2 = convert_bbox_to_geom_internal(bbox);

    initGEOS(lwnotice, lwgeom_geos_error);
    g1 = (GEOSGeometry *)LWGEOM2GEOS(geom1, LW_FALSE);
    g2 = (GEOSGeometry *)LWGEOM2GEOS(geom2, LW_FALSE);

    result = GEOSRelatePattern(g2, g1, "T*F*F*T**");
    GEOSGeom_destroy(g1);
    GEOSGeom_destroy(g2);
    lwgeom_free(geom1);
    lwgeom_free(geom2);
    return result;
}

/*is n_corner covered by bbox?*/
bool coveredBy(const N_CORNER *n_corner, const BBox *bbox)
{
    int result;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    LWGEOM *geom1 = convert_n_corner_to_geom_internal(n_corner);
    LWGEOM *geom2 = convert_bbox_to_geom_internal(bbox);

    initGEOS(lwnotice, lwgeom_geos_error);
    g1 = (GEOSGeometry *)LWGEOM2GEOS(geom1, LW_FALSE);
    g2 = (GEOSGeometry *)LWGEOM2GEOS(geom2, LW_FALSE);

    result = GEOSRelatePattern(g1, g2, "T*F*T*T**");
    GEOSGeom_destroy(g1);
    GEOSGeom_destroy(g2);
    lwgeom_free(geom1);
    lwgeom_free(geom2);
    return result;
}

/*does n_corner cover bbox?*/
bool covers(const N_CORNER *n_corner, const BBox *bbox)
{
    int result;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    LWGEOM *geom1 = convert_n_corner_to_geom_internal(n_corner);
    LWGEOM *geom2 = convert_bbox_to_geom_internal(bbox);

    initGEOS(lwnotice, lwgeom_geos_error);
    g1 = (GEOSGeometry *)LWGEOM2GEOS(geom1, LW_FALSE);
    g2 = (GEOSGeometry *)LWGEOM2GEOS(geom2, LW_FALSE);

    result = GEOSRelatePattern(g2, g1, "T*F*T*T**");
    GEOSGeom_destroy(g1);
    GEOSGeom_destroy(g2);
    lwgeom_free(geom1);
    lwgeom_free(geom2);
    return result;
}

/* this overlap is according to the 9-IM! */
bool overlap(const N_CORNER *n_corner, const BBox *bbox)
{
    int result;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    LWGEOM *geom1 = convert_n_corner_to_geom_internal(n_corner);
    LWGEOM *geom2 = convert_bbox_to_geom_internal(bbox);

    initGEOS(lwnotice, lwgeom_geos_error);
    g1 = (GEOSGeometry *)LWGEOM2GEOS(geom1, LW_FALSE);
    g2 = (GEOSGeometry *)LWGEOM2GEOS(geom2, LW_FALSE);

    result = GEOSOverlaps(g1, g2);
    GEOSGeom_destroy(g1);
    GEOSGeom_destroy(g2);
    lwgeom_free(geom1);
    lwgeom_free(geom2);
    return result;
}

/* this meet is according to the 9-IM! */
bool meet(const N_CORNER *n_corner, const BBox *bbox)
{
    int result;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    LWGEOM *geom1 = convert_n_corner_to_geom_internal(n_corner);
    LWGEOM *geom2 = convert_bbox_to_geom_internal(bbox);

    initGEOS(lwnotice, lwgeom_geos_error);
    g1 = (GEOSGeometry *)LWGEOM2GEOS(geom1, LW_FALSE);
    g2 = (GEOSGeometry *)LWGEOM2GEOS(geom2, LW_FALSE);

    result = GEOSTouches(g1, g2);
    GEOSGeom_destroy(g1);
    GEOSGeom_destroy(g2);
    lwgeom_free(geom1);
    lwgeom_free(geom2);
    return result;
}

bool equal(const N_CORNER *n_corner, const BBox *bbox)
{
    int result;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    LWGEOM *geom1;
    LWGEOM *geom2;

    if (n_corner->num > 4)
    {
        return false;
    }

    geom1 = convert_n_corner_to_geom_internal(n_corner);
    geom2 = convert_bbox_to_geom_internal(bbox);

    initGEOS(lwnotice, lwgeom_geos_error);
    g1 = (GEOSGeometry *)LWGEOM2GEOS(geom1, LW_FALSE);
    g2 = (GEOSGeometry *)LWGEOM2GEOS(geom2, LW_FALSE);

    result = GEOSEquals(g1, g2);
    GEOSGeom_destroy(g1);
    GEOSGeom_destroy(g2);
    lwgeom_free(geom1);
    lwgeom_free(geom2);
    return result;
}

bool contain_or_covers(const N_CORNER *n_corner, const BBox *bbox)
{
    int result;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    LWGEOM *geom1 = convert_n_corner_to_geom_internal(n_corner);
    LWGEOM *geom2 = convert_bbox_to_geom_internal(bbox);

    initGEOS(lwnotice, lwgeom_geos_error);
    g1 = (GEOSGeometry *)LWGEOM2GEOS(geom1, LW_FALSE);
    g2 = (GEOSGeometry *)LWGEOM2GEOS(geom2, LW_FALSE);

    result = GEOSContains(g1, g2);
    GEOSGeom_destroy(g1);
    GEOSGeom_destroy(g2);
    lwgeom_free(geom1);
    lwgeom_free(geom2);
    return result && !equal(n_corner, bbox);
}

bool n_corner_check_predicate(const SpatialApproximation *ap, const SpatialApproximation *bbox, uint8_t predicate)
{
    if (bbox->type == BBOX_TYPE)
    {
        N_CORNER_APPROX *n_corner_ap = (void *)ap;
        BBOX_APPROX *bbox_ap = (void *)bbox;
        N_CORNER *n_corner = &(n_corner_ap->n_corner);
        BBox *bbox1 = &(bbox_ap->bbox);

        switch (predicate)
        {
        case INTERSECTS:
            return intersect(n_corner, bbox1);
        case DISJOINT:
            return !intersect(n_corner, bbox1);
        case OVERLAP:
            return overlap(n_corner, bbox1);
        case MEET:
            return meet(n_corner, bbox1);
        case INSIDE:
            return inside(n_corner, bbox1);
        case CONTAINS:
            return contains(n_corner, bbox1);
        case COVEREDBY:
            return coveredBy(n_corner, bbox1);
        case COVERS:
            return covers(n_corner, bbox1);
        case EQUAL:
            return equal(n_corner, bbox1);
        case INSIDE_OR_COVEREDBY:
            return inside_or_coveredBy(n_corner, bbox1);
        case CONTAINS_OR_COVERS:
            return contain_or_covers(n_corner, bbox1);
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

size_t serialize_n_corner(const SpatialApproximation *ap, uint8_t *buf)
{
    uint8_t *loc;
    uint8_t n;
    int i, max = 5;
    N_CORNER_APPROX *n_corner_ap = (void *)ap;
    ApproxPoint *p;
    ApproxPoint nil = {0, 0};
    loc = buf;

    memcpy(loc, &(ap->id), sizeof(uint64_t));
    loc += sizeof(uint64_t);

    memcpy(loc, &(ap->type), sizeof(uint8_t));
    loc += sizeof(uint8_t);

    p = n_corner_ap->n_corner.p;
    n = n_corner_ap->n_corner.num;

    memcpy(loc, &n, sizeof(uint8_t));
    loc += sizeof(uint8_t);

    for (i = 0; i < n; i++)
    {
        memcpy(loc, &(p[i]), sizeof(ApproxPoint));
        loc += sizeof(ApproxPoint);
    }

    if (ap->type == N_CORNER_4_TYPE)
        max = 4;

    for (i = 0; i < (max - n); i++)
    {
        memcpy(loc, &nil, sizeof(ApproxPoint));
        loc += sizeof(ApproxPoint);
    }

    return (size_t)(loc - buf);
}

size_t n_5_corner_size()
{
    return (sizeof(ApproxPoint) * 5) + sizeof(uint64_t) + (sizeof(uint8_t) * 2);
}

size_t n_4_corner_size()
{
    return (sizeof(ApproxPoint) * 4) + sizeof(uint64_t) + (sizeof(uint8_t) * 2);
}

size_t n_corner_size(uint8_t type)
{
    if (type == N_CORNER_5_TYPE)
        return n_5_corner_size();
    else
        return n_4_corner_size();
}

bool write_n_corner_to_file(const FileSpecification *fs, const SpatialApproximation **ap, int n)
{
    size_t total_size = (n * n_corner_size(ap[0]->type)) + sizeof(int), finalsize;
    uint8_t *buf, *loc;
    int i;

    if (total_size > fs->page_size)
    {
        _DEBUG(ERROR, "Total size is bigger than page size\n");
        return false;
    }

    if (fs->io_access == DIRECT_ACCESS)
    {
        if (posix_memalign((void **)&buf, fs->page_size, fs->page_size))
        {
            _DEBUG(ERROR, "Allocation failed at N_CORNER write_to_file\n");
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
        loc += serialize_n_corner(ap[i], loc);
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

SpatialApproximation *recreate_n_corner(int n, ApproxPoint *p, uint8_t type, uint64_t id)
{
    static const SpatialApproximationInterface vtable = {convert_n_corner_to_geom, n_corner_check_predicate, destroy_n_corner};
    static SpatialApproximation base = {&vtable};
    N_CORNER_APPROX *approx;
    int i;

    N_CORNER n_corner;
    n_corner.num = n;
    n_corner.p = n_corner_allocate_points(n);

    for (i = 0; i < n; i++)
    {
        n_corner.p[i] = p[i];
    }

    base.type = type;
    base.id = id;
    approx = lwalloc(sizeof(N_CORNER_APPROX));
    memcpy(&approx->base, &base, sizeof(base));
    approx->n_corner = n_corner;

    return &approx->base;
}

SpatialApproximation *read_serialized_n_corner(uint8_t *buf, size_t *size)
{
    SpatialApproximation *ap;
    uint8_t *start_ptr = buf;
    uint8_t n, type;
    uint64_t id;
    ApproxPoint p[5];
    int i, max = 5;

    memcpy(&id, buf, sizeof(uint64_t));
    buf += sizeof(uint64_t);

    memcpy(&type, buf, sizeof(uint8_t));
    buf += sizeof(uint8_t);

    memcpy(&n, buf, sizeof(uint8_t));
    buf += sizeof(uint8_t);

    for (i = 0; i < n; i++)
    {
        memcpy(&(p[i]), buf, sizeof(ApproxPoint));
        buf += sizeof(ApproxPoint);
    }

    if (type == N_CORNER_4_TYPE)
        max = 4;

    buf += sizeof(ApproxPoint) * (max - n);

    ap = recreate_n_corner(n, p, type, id);

    *size = (size_t)(buf - start_ptr);
    return ap;
}

SpatialApproximation **read_n_corner_from_file(const FileSpecification *fs, int page, int *n)
{
    uint8_t *buf, *loc;
    SpatialApproximation **ap;
    int i;
    size_t size;

    if (fs->io_access == DIRECT_ACCESS)
    {
        if (posix_memalign((void **)&buf, fs->page_size, fs->page_size))
        {
            _DEBUG(ERROR, "Allocation failed at N_CORNER read_from_file\n");
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
        ap[i] = read_serialized_n_corner(loc, &size);
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

SpatialApproximation *create_n_corner(LWGEOM *geom, uint8_t type, uint64_t id)
{
    static const SpatialApproximationInterface vtable = {convert_n_corner_to_geom, n_corner_check_predicate, destroy_n_corner};
    static SpatialApproximation base = {&vtable};
    N_CORNER_APPROX *approx;

    N_CORNER n_corner;
    if (type == N_CORNER_4_TYPE)
    {
        convert_geom_to_4_corner(geom, &n_corner);
    }
    else if (type == N_CORNER_5_TYPE)
    {
        convert_geom_to_5_corner(geom, &n_corner);
    }

    base.type = type;
    base.id = id;
    approx = lwalloc(sizeof(N_CORNER_APPROX));
    memcpy(&approx->base, &base, sizeof(base));
    approx->n_corner = n_corner;

    return &approx->base;
}

void destroy_n_corner(SpatialApproximation *ap)
{
    N_CORNER_APPROX *n_corner_ap = (void *)ap;
    n_corner_free(&n_corner_ap->n_corner);
    lwfree(n_corner_ap);
}