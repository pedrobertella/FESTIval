#ifndef _RMBR_HANDLER_H
#define _RMBR_HANDLER_H

#include <liblwgeom.h>
#include <stdbool.h>
#include "SpatialApproximation.h"
#include "bbox_approx_handler.h"
#include "RMBR.h"

typedef struct {
    SpatialApproximation base; 
    RMBR rmbr;
} RMBR_APPROX;

extern SpatialApproximation *create_rmbr(LWGEOM *geom, uint64_t id);

extern size_t rmbr_size(void);

extern bool write_rmbr_to_file(const FileSpecification *fs, const SpatialApproximation **ap, int n);

extern SpatialApproximation *read_serialized_rmbr(uint8_t *buf, size_t *size);

extern SpatialApproximation **read_rmbr_from_file(const FileSpecification *fs, int page, int *n);

#endif /* _RMBR_HANDLER_H */