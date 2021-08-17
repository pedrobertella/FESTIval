#ifndef _RMBP_HANDLER_H
#define _RMBP_HANDLER_H

#include <liblwgeom.h>
#include <stdbool.h>
#include "SpatialApproximation.h"
#include "bbox_approx_handler.h"
#include "RMBP.h"

typedef struct {
    SpatialApproximation base; 
    RMBP rmbp;
} RMBP_APPROX;

extern SpatialApproximation *create_rmbp(LWGEOM *geom, uint64_t id);

extern size_t rmbp_size(void);

extern bool write_rmbp_to_file(const FileSpecification *fs, const SpatialApproximation **ap, int n);

extern SpatialApproximation *read_serialized_rmbp(uint8_t *buf, size_t *size);

extern SpatialApproximation **read_rmbp_from_file(const FileSpecification *fs, int page, int *n);

#endif /* _RMBP_HANDLER_H */