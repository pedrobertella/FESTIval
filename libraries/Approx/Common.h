#ifndef _COMMON_H
#define _COMMON_H
#include <stdio.h>

typedef struct point2d {
	double x, y;
} ApproxPoint;

typedef enum {
	EXTERIOR = -1,
	BOUNDARY,
	INTERIOR
} POLY_REGION;

#endif