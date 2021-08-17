#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif
#include "Common.h"

typedef void* RMBP;

EXTERNC RMBP rmbp_create(ApproxPoint a, ApproxPoint b, ApproxPoint c, ApproxPoint d);
EXTERNC void rmbp_free(RMBP rmbp);
EXTERNC bool rmbp_is_in_exterior(RMBP rmbp, ApproxPoint p);
EXTERNC bool rmbp_is_in_boundary(RMBP rmbp, ApproxPoint p);
EXTERNC bool rmbp_is_in_interior(RMBP rmbp, ApproxPoint p);
EXTERNC int rmbp_region(RMBP rmbp, ApproxPoint p);
EXTERNC RMBP rmbp_calculate(int n, ApproxPoint* P);
EXTERNC ApproxPoint *rmbp_get_points(RMBP rmbp);