/**********************************************************************
 *
 * FESTIval - Framework to Evaluate SpaTial Indices in non-VolAtiLe memories and hard disk drives.
 * http://gbd.dc.ufscar.br/festival/
 *
 * Copyright (C) 2016-2018 Anderson Chaves Carniel <accarniel@gmail.com>
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU General Public Licence. See the COPYING file.
 *
 * Fully developed by Anderson Chaves Carniel
 *
 **********************************************************************/

#include <string.h>
#include <stdbool.h>
#include <liblwgeom.h>
#include "bbox_approx_handler.h"
#include "../main/math_util.h"
#include "../main/log_messages.h"

static void convert_geom_to_bbox(const LWGEOM *geom, BBox *bbox);
static LWGEOM *convert_bbox_to_geom(const SpatialApproximation *ap);
static bool bbox_approx_check_predicate(const SpatialApproximation *ap, const SpatialApproximation *bbox, uint8_t predicate);
static size_t serialize_bbox(const SpatialApproximation *ap, uint8_t *buf);
static void destroy_bbox(SpatialApproximation *ap);

void convert_geom_to_bbox(const LWGEOM *geom, BBox *bbox)
{
    GBOX gbox;
    lwgeom_calculate_gbox(geom, &gbox);

    if (NUM_OF_DIM == 2)
    {
        bbox->min[0] = gbox.xmin;
        bbox->min[1] = gbox.ymin;

        bbox->max[0] = gbox.xmax;
        bbox->max[1] = gbox.ymax;
    }
}

LWGEOM *convert_bbox_to_geom(const SpatialApproximation *ap)
{
    BBOX_APPROX *approx = (void *)ap;
    BBox *box = &(approx->bbox);

    if (NUM_OF_DIM == 2)
    {
        /*
         * In order to always return a valid geometry:
         *     - If the bounding box is a single point then return a
         *     POINT geometry
         *     - If the bounding box represents either a horizontal or
         *     vertical line, return a LINESTRING geometry
         *     - Otherwise return a POLYGON
         */

        if ((box->min[0] == box->max[0]) && (box->min[1] == box->max[1]))
        {
            /* Construct a point */
            LWPOINT *point = lwpoint_make2d(SRID_UNKNOWN, box->min[0], box->min[1]);
            return lwpoint_as_lwgeom(point);
        }
        else if ((box->min[0] == box->max[0]) || (box->min[1] == box->max[1]))
        {
            LWLINE *line;
            POINTARRAY *pa;
            POINT4D pt;

            /* Construct point array */
            pa = ptarray_construct_empty(0, 0, 2);

            /* Assign coordinates to POINT2D array */
            pt.x = box->min[0]; //xmin;
            pt.y = box->min[1]; //ymin;
            ptarray_append_point(pa, &pt, LW_TRUE);
            pt.x = box->max[0]; //xmax;
            pt.y = box->max[1]; //ymax;
            ptarray_append_point(pa, &pt, LW_TRUE);

            /* Construct and serialize linestring */
            line = lwline_construct(SRID_UNKNOWN, NULL, pa);
            return lwline_as_lwgeom(line);
        }
        else
        {
            LWPOLY *poly;
            POINTARRAY *pa = ptarray_construct_empty(0, 0, 5);
            POINT4D pt;
            POINTARRAY **ppa = lwalloc(sizeof(POINTARRAY *));

            /* Assign coordinates to point array */
            pt.x = box->min[0]; //xmin
            pt.y = box->min[1]; //ymin
            ptarray_append_point(pa, &pt, LW_TRUE);
            pt.x = box->min[0]; //xmin
            pt.y = box->max[1]; //ymax
            ptarray_append_point(pa, &pt, LW_TRUE);
            pt.x = box->max[0]; //xmax
            pt.y = box->max[1]; //ymax
            ptarray_append_point(pa, &pt, LW_TRUE);
            pt.x = box->max[0]; //xmax
            pt.y = box->min[1]; //ymin
            ptarray_append_point(pa, &pt, LW_TRUE);
            pt.x = box->min[0]; //xmin
            pt.y = box->min[1]; //ymin
            ptarray_append_point(pa, &pt, LW_TRUE);

            /* Construct polygon */
            ppa[0] = pa;
            poly = lwpoly_construct(SRID_UNKNOWN, NULL, 1, ppa);
            return lwpoly_as_lwgeom(poly);
        }
    }
    else
    {
        return NULL;
    }
}

LWGEOM *convert_bbox_to_geom_internal(const BBox *box)
{
    if (NUM_OF_DIM == 2)
    {
        /*
         * In order to always return a valid geometry:
         *     - If the bounding box is a single point then return a
         *     POINT geometry
         *     - If the bounding box represents either a horizontal or
         *     vertical line, return a LINESTRING geometry
         *     - Otherwise return a POLYGON
         */

        if ((box->min[0] == box->max[0]) && (box->min[1] == box->max[1]))
        {
            /* Construct a point */
            LWPOINT *point = lwpoint_make2d(SRID_UNKNOWN, box->min[0], box->min[1]);
            return lwpoint_as_lwgeom(point);
        }
        else if ((box->min[0] == box->max[0]) || (box->min[1] == box->max[1]))
        {
            LWLINE *line;
            POINTARRAY *pa;
            POINT4D pt;

            /* Construct point array */
            pa = ptarray_construct_empty(0, 0, 2);

            /* Assign coordinates to POINT2D array */
            pt.x = box->min[0]; //xmin;
            pt.y = box->min[1]; //ymin;
            ptarray_append_point(pa, &pt, LW_TRUE);
            pt.x = box->max[0]; //xmax;
            pt.y = box->max[1]; //ymax;
            ptarray_append_point(pa, &pt, LW_TRUE);

            /* Construct and serialize linestring */
            line = lwline_construct(SRID_UNKNOWN, NULL, pa);
            return lwline_as_lwgeom(line);
        }
        else
        {
            LWPOLY *poly;
            POINTARRAY *pa = ptarray_construct_empty(0, 0, 5);
            POINT4D pt;
            POINTARRAY **ppa = lwalloc(sizeof(POINTARRAY *));

            /* Assign coordinates to point array */
            pt.x = box->min[0]; //xmin
            pt.y = box->min[1]; //ymin
            ptarray_append_point(pa, &pt, LW_TRUE);
            pt.x = box->min[0]; //xmin
            pt.y = box->max[1]; //ymax
            ptarray_append_point(pa, &pt, LW_TRUE);
            pt.x = box->max[0]; //xmax
            pt.y = box->max[1]; //ymax
            ptarray_append_point(pa, &pt, LW_TRUE);
            pt.x = box->max[0]; //xmax
            pt.y = box->min[1]; //ymin
            ptarray_append_point(pa, &pt, LW_TRUE);
            pt.x = box->min[0]; //xmin
            pt.y = box->min[1]; //ymin
            ptarray_append_point(pa, &pt, LW_TRUE);

            /* Construct polygon */
            ppa[0] = pa;
            poly = lwpoly_construct(SRID_UNKNOWN, NULL, 1, ppa);
            return lwpoly_as_lwgeom(poly);
        }
    }
    else
    {
        return NULL;
    }
}

bool bbox_approx_check_predicate(const SpatialApproximation *ap, const SpatialApproximation *bbox, uint8_t predicate)
{
    if (bbox->type == BBOX_TYPE)
    {
        BBOX_APPROX *bbox1_ap = (void *)ap;
        BBOX_APPROX *bbox2_ap = (void *)bbox;
        BBox *bbox1 = &(bbox1_ap->bbox);
        BBox *bbox2 = &(bbox2_ap->bbox);

        return bbox_check_predicate(bbox1, bbox2, predicate);
    }
    else
    {
        _DEBUG(WARNING, "Wrong type for bbox approximation\n");
        return false;
    }
}

size_t serialize_bbox(const SpatialApproximation *ap, uint8_t *buf)
{
    uint8_t *loc;
    int i;
    BBOX_APPROX *bbox_ap;
    loc = buf;

    bbox_ap = (void *)ap;

    memcpy(loc, &(ap->id), sizeof(uint64_t));
    loc += sizeof(uint64_t);

    for (i = 0; i < NUM_OF_DIM; i++)
    {
        memcpy(loc, &(bbox_ap->bbox.min[i]), sizeof(double));
        loc += sizeof(double);
    }

    for (i = 0; i < NUM_OF_DIM; i++)
    {
        memcpy(loc, &(bbox_ap->bbox.max[i]), sizeof(double));
        loc += sizeof(double);
    }

    return (size_t)(loc - buf);
}

size_t bbox_size()
{
    return (sizeof(double) * NUM_OF_DIM * 2) + sizeof(uint64_t);
}

bool write_bbox_to_file(const FileSpecification *fs, const SpatialApproximation **ap, int n)
{
    uint8_t *buf, *loc;
    size_t total_size = (n * bbox_size()) + sizeof(int), finalsize;
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
            _DEBUG(ERROR, "Allocation failed at BBOX write_to_file\n");
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
        loc += serialize_bbox(ap[i], loc);
    }

    finalsize = (size_t)(loc - buf);

    if (total_size != finalsize) /* Uh oh! */
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

SpatialApproximation *recreate_bbox(double *min, double *max, uint64_t id)
{
    static const SpatialApproximationInterface vtable = {convert_bbox_to_geom, bbox_approx_check_predicate, destroy_bbox};
    static SpatialApproximation base = {&vtable};
    BBOX_APPROX *approx;

    BBox bbox;
    int i;
    for (i = 0; i < NUM_OF_DIM; i++)
    {
        bbox.min[i] = min[i];
        bbox.max[i] = max[i];
    }

    base.type = BBOX_TYPE;
    base.id = id;
    approx = lwalloc(sizeof(BBOX_APPROX));
    memcpy(&approx->base, &base, sizeof(base));
    approx->bbox = bbox;

    return &approx->base;
}

SpatialApproximation *read_serialized_bbox(uint8_t *buf, size_t *size)
{
    uint8_t *start_ptr = buf;
    double min[NUM_OF_DIM], max[NUM_OF_DIM];
    uint64_t id;
    int i;
    SpatialApproximation *ap;

    memcpy(&id, buf, sizeof(uint64_t));
    buf += sizeof(uint64_t);

    for (i = 0; i < NUM_OF_DIM; i++)
    {
        memcpy(&(min[i]), buf, sizeof(double));
        buf += sizeof(double);
    }

    for (i = 0; i < NUM_OF_DIM; i++)
    {
        memcpy(&(max[i]), buf, sizeof(double));
        buf += sizeof(double);
    }

    ap = recreate_bbox(min, max, id);

    *size = (size_t)(buf - start_ptr);
    return ap;
}

SpatialApproximation **read_bbox_from_file(const FileSpecification *fs, int page, int *n)
{
    uint8_t *buf, *loc;
    SpatialApproximation **ap;
    int i;
    size_t size;

    if (fs->io_access == DIRECT_ACCESS)
    {
        if (posix_memalign((void **)&buf, fs->page_size, fs->page_size))
        {
            _DEBUG(ERROR, "Allocation failed at MBC read_from_file\n");
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
        ap[i] = read_serialized_bbox(loc, &size);
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

SpatialApproximation *create_bbox_approx(LWGEOM *geom, uint64_t id)
{
    static const SpatialApproximationInterface vtable = {convert_bbox_to_geom, bbox_approx_check_predicate, destroy_bbox};
    static SpatialApproximation base = {&vtable};
    BBOX_APPROX *approx;
    
    BBox bbox;
    convert_geom_to_bbox(geom, &bbox);
    base.type = BBOX_TYPE;
    base.id = id;
    approx = lwalloc(sizeof(BBOX_APPROX));
    memcpy(&approx->base, &base, sizeof(base));
    approx->bbox = bbox;

    return &approx->base;
}

void destroy_bbox(SpatialApproximation *ap)
{
    BBOX_APPROX *bbox_ap = (void *)ap;
    lwfree(bbox_ap);
}