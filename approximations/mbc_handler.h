#ifndef _MBC_HANDLER_H
#define _MBC_HANDLER_H

#include <liblwgeom.h>
#include <stdbool.h>
#include "SpatialApproximation.h"
#include "bbox_approx_handler.h"

typedef struct Mbc {
	double c[NUM_OF_DIM];
	double r;
} MBC;

typedef struct {
    SpatialApproximation base; 
    MBC mbc;
} MBC_APPROX;

extern SpatialApproximation *create_mbc(LWGEOM *geom, uint64_t id);

extern size_t mbc_size(void);

extern bool write_mbc_to_file(const FileSpecification *fs, const SpatialApproximation **ap, int n);

extern SpatialApproximation *read_serialized_mbc(uint8_t *buf, size_t *size);

extern SpatialApproximation **read_mbc_from_file(const FileSpecification *fs, int page, int *n);

#endif /* _MBC_HANDLER_H */
