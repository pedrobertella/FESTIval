#ifndef _N_CORNER_HANDLER_H
#define _N_CORNER_HANDLER_H

#include <liblwgeom.h>
#include <stdbool.h>
#include "SpatialApproximation.h"
#include "bbox_approx_handler.h"
#include "N_CORNER.h"

typedef struct {
    SpatialApproximation base; 
    N_CORNER n_corner;
} N_CORNER_APPROX;

extern SpatialApproximation *create_n_corner(LWGEOM *geom, uint8_t type, uint64_t id);

extern size_t n_5_corner_size(void);

extern size_t n_4_corner_size(void);

extern bool write_n_corner_to_file(const FileSpecification *fs, const SpatialApproximation **ap, int n);

extern SpatialApproximation *read_serialized_n_corner(uint8_t *buf, size_t *size);

extern SpatialApproximation **read_n_corner_from_file(const FileSpecification *fs, int page, int *n);

#endif /* _N_CORNER_HANDLER_H */