#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif
#include "Common.h"

typedef void* RMBR;

EXTERNC RMBR rmbr_create(ApproxPoint a, ApproxPoint b, ApproxPoint c, ApproxPoint d);
EXTERNC void rmbr_free(RMBR rmbr);
EXTERNC bool rmbr_is_in_exterior(RMBR rmbr, ApproxPoint p);
EXTERNC bool rmbr_is_in_boundary(RMBR rmbr, ApproxPoint p);
EXTERNC bool rmbr_is_in_interior(RMBR rmbr, ApproxPoint p);
EXTERNC int rmbr_region(RMBR rmbr, ApproxPoint p);
EXTERNC RMBR rmbr_calculate(int n, ApproxPoint *P);
EXTERNC ApproxPoint *rmbr_get_points(RMBR rmbr);