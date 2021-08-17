#ifndef _LWGEOM_HANDLER_H
#define _LWGEOM_HANDLER_H

#include <liblwgeom.h>
#include <stdbool.h>
#include <geos_c.h>

extern void lwgeom_force_counterclockwise(LWGEOM *geom);

extern LWGEOM *convexhull(LWGEOM *geom);

extern void lwgeom_free_array(int n, LWGEOM **geoms);

#endif /* _LWGEOM_HANDLER_H */