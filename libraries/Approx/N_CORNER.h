#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

#include "Common.h"

typedef struct {
    ApproxPoint *p;
    int num;
} N_CORNER;

EXTERNC N_CORNER five_corner_calculate(int n, ApproxPoint* P);
EXTERNC N_CORNER four_corner_calculate(int n, ApproxPoint* P);

EXTERNC ApproxPoint* n_corner_allocate_points(int n);
EXTERNC void n_corner_free(N_CORNER *n_corner);