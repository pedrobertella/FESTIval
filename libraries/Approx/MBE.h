#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

#include "Common.h"

typedef void* MBE;

EXTERNC MBE mbe_create(int n, ApproxPoint *P);
EXTERNC void mbe_free(MBE mbe);
EXTERNC bool mbe_is_in_exterior(MBE mbe, ApproxPoint p);
EXTERNC bool mbe_is_in_boundary(MBE mbe, ApproxPoint p);
EXTERNC bool mbe_is_in_interior(MBE mbe, ApproxPoint p);
EXTERNC int mbe_region(MBE mbe, ApproxPoint p);
EXTERNC MBE mbe_calculate(int n, ApproxPoint* P);
EXTERNC ApproxPoint *mbe_get_points(int *n, MBE mbe);