#ifndef _SPATIAL_APPROXIMATION_H
#define _SPATIAL_APPROXIMATION_H

#include <stdio.h>
#include <stdbool.h>
#include <malloc.h> 
#include <liblwgeom.h>
#include "../main/io_handler.h"

#define NUM_OF_DIM 2
#define MAX_DIM ((NUM_OF_DIM)-1)

/* Predicates */

#define INTERSECTS 1
#define OVERLAP 2
#define DISJOINT 3
#define MEET 4
#define INSIDE 5
#define COVEREDBY 6
#define CONTAINS 7
#define COVERS 8
#define EQUAL 9
#define INSIDE_OR_COVEREDBY 10
#define CONTAINS_OR_COVERS 11

/* Spatial approximation types */

#define BBOX_TYPE 1
#define MBC_TYPE 2
#define RMBR_TYPE 3
#define RMBP_TYPE 4
#define MBE_TYPE 5
#define N_CORNER_TYPE 6
#define N_CORNER_4_TYPE 0b00000110
#define N_CORNER_5_TYPE 0b00010110

/* Struct definitions */

typedef struct
{
    const struct _SpatialApproximationInterface *const vtable;
    uint8_t type;
    uint64_t id;
} SpatialApproximation;

typedef struct _SpatialApproximationInterface
{
    LWGEOM *(*convert_to_lwgeom)(const SpatialApproximation *ap);
    bool (*check_predicate)(const SpatialApproximation *ap, const SpatialApproximation *bbox, uint8_t predicate);
    void (*free_approx)(SpatialApproximation *ap);
} SpatialApproximationInterface;

/* Wrapper functions */

static inline LWGEOM *spatialapproximation_convert_to_lwgeom(const SpatialApproximation *ap)
{
    return ap->vtable->convert_to_lwgeom(ap);
}

static inline bool spatialapproximation_check_predicate(const SpatialApproximation *ap, const SpatialApproximation *bbox, uint8_t p)
{
    return ap->vtable->check_predicate(ap, bbox, p);
}

static inline void spatialapproximation_free(SpatialApproximation *ap){
    ap->vtable->free_approx(ap);
}

#endif /* _SPATIAL_APPROXIMATION_H */