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

/* 
 * File:   bbox_handler.h
 * Author: Anderson Chaves Carniel
 *
 * Created on February 25, 2016, 7:57 PM
 */

#ifndef _BBOX_APPROX_HANDLER_H
#define _BBOX_APPROX_HANDLER_H

#include <liblwgeom.h>
#include <stdbool.h>
#include "SpatialApproximation.h"
#include "../main/bbox_handler.h"

typedef struct {
    SpatialApproximation base; 
    BBox bbox;
} BBOX_APPROX;

extern SpatialApproximation *create_bbox_approx(LWGEOM *geom, uint64_t id);

extern LWGEOM *convert_bbox_to_geom_internal(const BBox *box);

extern size_t bbox_size(void);

extern bool write_bbox_to_file(const FileSpecification *fs, const SpatialApproximation **ap, int n);

extern SpatialApproximation *read_serialized_bbox(uint8_t *buf, size_t *size);

extern SpatialApproximation **read_bbox_from_file(const FileSpecification *fs, int page, int *n);

extern SpatialApproximation *recreate_bbox(double *min, double *max, uint64_t id);

#endif /* _BBOX_APPROX_HANDLER_H */

