#ifndef _MBE_HANDLER_H
#define _MBE_HANDLER_H

#include <liblwgeom.h>
#include <stdbool.h>
#include "SpatialApproximation.h"
#include "bbox_approx_handler.h"
#include "MBE.h"

typedef struct {
    SpatialApproximation base; 
    MBE mbe;
} MBE_APPROX;

extern SpatialApproximation *create_mbe(LWGEOM *geom, uint64_t id);

extern size_t mbe_size(void);

extern bool write_mbe_to_file(const FileSpecification *fs, const SpatialApproximation **ap, int n);

extern SpatialApproximation *read_serialized_mbe(uint8_t *buf, size_t *size);

extern SpatialApproximation **read_mbe_from_file(const FileSpecification *fs, int page, int *n);

#endif /* _MBE_HANDLER_H */