#include "lwgeom_handler.h"
#include "lwgeom_log.h"
#include "lwgeom_geos.h"

extern void lwgeom_force_counterclockwise(LWGEOM *geom)
{
    lwgeom_force_clockwise(geom);
    lwgeom_reverse_in_place(geom);
}

extern LWGEOM *convexhull(LWGEOM *geom)
{
    LWGEOM *result;
    GEOSGeometry *g1, *g3;
    initGEOS(lwnotice, lwgeom_geos_error);

    g1 = (GEOSGeometry *)LWGEOM2GEOS(geom, LW_FALSE);

    if (0 == g1)
    { /* exception thrown at construction */
        lwnotice("First argument geometry could not be converted to GEOS: %s", lwgeom_geos_errmsg);
        return NULL;
    }

    g3 = (GEOSGeometry *)GEOSConvexHull(g1);
    GEOSGeom_destroy(g1);

    if (g3 == NULL)
    {
        lwnotice("GEOSConvexHull: %s", lwgeom_geos_errmsg);
        return NULL; /* never get here */
    }

    GEOSSetSRID(g3, geom->srid);

    result = GEOS2LWGEOM(g3, 0);
    GEOSGeom_destroy(g3);

    if (result == NULL)
    {
        return NULL;
    }

    return result;
}

void lwgeom_free_array(int n, LWGEOM **geoms) {
	int i;
	for (i = 0; i < n; ++i)
		lwgeom_free (geoms[i]);
	lwfree (geoms);
}