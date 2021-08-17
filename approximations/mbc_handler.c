#include <stdlib.h>
#include <stdio.h>
#include <liblwgeom.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "../main/math_util.h"
#include "mbc_handler.h"
#include "../main/log_messages.h"

static double euclidean_distance(const double *a, const double *b);
static bool intersect_from_bbox(const MBC *mbc, const BBox *bbox);
static bool inside_from_bbox(const MBC *mbc, const BBox *bbox);
static bool contains_from_bbox(const MBC *mbc, const BBox *bbox);
static bool coveredBy_from_bbox(const MBC *mbc, const BBox *bbox);
static bool covers_from_bbox(const MBC *mbc, const BBox *bbox);
static bool equal_from_bbox(const MBC *mbc, const BBox *bbox);
static bool inside_or_coveredBy_from_bbox(const MBC *mbc, const BBox *bbox);
static bool contain_or_covers_from_bbox(const MBC *mbc, const BBox *bbox);
static bool meet_from_bbox(const MBC *mbc, const BBox *bbox);
static bool overlap_from_bbox(const MBC *mbc, const BBox *bbox);
static void convert_geom_to_mbc(const LWGEOM *lw, MBC *mbc);
static LWGEOM *convert_mbc_to_geom(const SpatialApproximation *ap);
static bool mbc_check_predicate_from_bbox(const SpatialApproximation *ap, const SpatialApproximation *bbox, uint8_t predicate);
static size_t serialize_mbc(const SpatialApproximation *ap, uint8_t *buf);
static SpatialApproximation *recreate_mbc(double r, double *c, uint64_t id);
static void destroy_mbc(SpatialApproximation *ap);

double euclidean_distance(const double *a, const double *b)
{
	return sqrt(pow(b[0] - a[0], 2) + pow(b[1] - a[1], 2));
}

void convert_geom_to_mbc(const LWGEOM *lw, MBC *mbc)
{
	LWBOUNDINGCIRCLE *aux = lwgeom_calculate_mbc(lw);

	if (NUM_OF_DIM == 2)
	{
		mbc->c[0] = aux->center->x;
		mbc->c[1] = aux->center->y;
		mbc->r = aux->radius;
	}

	lwboundingcircle_destroy(aux);
}

LWGEOM *convert_mbc_to_geom(const SpatialApproximation *ap)
{
	MBC_APPROX *approx = (void *)ap;
	MBC *mbc = &(approx->mbc);
	return lwpoly_as_lwgeom(lwpoly_construct_circle(SRID_DEFAULT, mbc->c[0], mbc->c[1], mbc->r, 50, LW_TRUE));
}

/* do mbc and bbox intersect? OK*/
bool intersect_from_bbox(const MBC *mbc, const BBox *bbox)
{
	double p1[2], p2[2];

	/*1- checking if mbc's center is inside bbox*/
	if ((DB_LE(bbox->min[0], mbc->c[0]) && DB_GE(bbox->max[0], mbc->c[0])) && (DB_LE(bbox->min[1], mbc->c[1]) && DB_GE(bbox->max[1], mbc->c[1])))
	{
		return true;
	}

	/*2- checking if each any point of mbc is intersecting bbox*/
	if (DB_LE(bbox->min[0], (mbc->c[0]) - mbc->r) && DB_GE(bbox->max[0], (mbc->c[0]) - mbc->r)) //left
	{
		if ((DB_LE(bbox->min[1], (mbc->c[1]) - mbc->r) && DB_GE(bbox->max[1], (mbc->c[1]) - mbc->r)) || (DB_LE(bbox->min[1], (mbc->c[1]) + mbc->r) && DB_GE(bbox->max[1], (mbc->c[1]) + mbc->r))) //bottom or top are in range
		{
			return true;
		}
	}
	if (DB_LE(bbox->min[0], (mbc->c[0]) + mbc->r) && DB_GE(bbox->max[0], (mbc->c[0]) + mbc->r)) //right
	{
		if ((DB_LE(bbox->min[1], (mbc->c[1]) - mbc->r) && DB_GE(bbox->max[1], (mbc->c[1]) - mbc->r)) || (DB_LE(bbox->min[1], (mbc->c[1]) + mbc->r) && DB_GE(bbox->max[1], (mbc->c[1]) + mbc->r))) //bottom or top are in range
		{
			return true;
		}
	}

	/*3- checking if each any of bbox point is intersecting mbc*/
	if (DB_LE(euclidean_distance(bbox->min, mbc->c), mbc->r))
	{
		return true;
	}
	if (DB_LE(euclidean_distance(bbox->max, mbc->c), mbc->r))
	{
		return true;
	}

	p1[0] = bbox->min[0];
	p1[1] = bbox->max[1];
	p2[0] = bbox->max[0];
	p2[1] = bbox->min[1];
	if (DB_LE(euclidean_distance(p1, mbc->c), mbc->r))
	{
		return true;
	}
	if (DB_LE(euclidean_distance(p2, mbc->c), mbc->r))
	{
		return true;
	}

	return false;
}

/* is mbc inside bbox? OK*/
bool inside_from_bbox(const MBC *mbc, const BBox *bbox)
{
	/*checking if each side of mbc is inside bbox*/
	if (DB_LT(bbox->min[0], (mbc->c[0]) - mbc->r) && DB_GT(bbox->max[0], (mbc->c[0]) + mbc->r))
	{
		if (DB_LT(bbox->min[1], (mbc->c[1]) - mbc->r) && DB_GT(bbox->max[1], (mbc->c[1]) + mbc->r))
		{
			return true;
		}
	}

	return false;
}

/* does mbc contain bbox? OK*/
bool contains_from_bbox(const MBC *mbc, const BBox *bbox)
{
	double p1[2], p2[2];

	/* checks if the distance between bbox points and mbc center is greater than or equal to mbcs radius, if so = false*/
	if (DB_GE(euclidean_distance(bbox->min, mbc->c), mbc->r))
	{
		return false;
	}
	if (DB_GE(euclidean_distance(bbox->max, mbc->c), mbc->r))
	{
		return false;
	}

	p1[0] = bbox->min[0];
	p1[1] = bbox->max[1];
	p2[0] = bbox->max[0];
	p2[1] = bbox->min[1];
	if (DB_GE(euclidean_distance(p1, mbc->c), mbc->r))
	{
		return false;
	}
	if (DB_GE(euclidean_distance(p2, mbc->c), mbc->r))
	{
		return false;
	}
	return true;
}

/* is mbc covered by bbox? */
bool coveredBy_from_bbox(const MBC *mbc, const BBox *bbox)
{
	/* checking for equal border points */
	int anyBorder = 0;
	if (DB_IS_EQUAL(bbox->min[0], (mbc->c[0]) - mbc->r) && (DB_GT(bbox->min[1], mbc->c[1]) && DB_LT(bbox->max[1], mbc->c[1]))) //left: mbc center minus radius X axis / y axis in inside bbox min and max
	{
		anyBorder++;
	}
	if (DB_IS_EQUAL(bbox->max[0], (mbc->c[0]) + mbc->r) && (DB_GT(bbox->min[1], mbc->c[1]) && DB_LT(bbox->max[1], mbc->c[1]))) //right: mbc center plus radius X axis / y axis in inside bbox min and max
	{
		anyBorder++;
	}

	if (DB_IS_EQUAL(bbox->min[1], (mbc->c[1]) - mbc->r) && (DB_GT(bbox->min[0], mbc->c[0]) && DB_LT(bbox->max[0], mbc->c[0]))) //bottom: mbc center minus radius Y axis / x axis in inside bbox min and max
	{
		anyBorder++;
	}
	if (DB_IS_EQUAL(bbox->max[1], (mbc->c[1]) + mbc->r) && (DB_GT(bbox->min[0], mbc->c[0]) && DB_LT(bbox->max[0], mbc->c[0]))) //bottom: mbc center plus radius Y axis / x axis in inside bbox min and max
	{
		anyBorder++;
	}
	if (anyBorder == 0)
	{
		return false;
	}
	/*checking if each side of mbc is inside bbox or on its border*/
	if (DB_LE(bbox->min[0], (mbc->c[0]) - mbc->r) && DB_GE(bbox->max[0], (mbc->c[0]) + mbc->r))
	{
		if (DB_LE(bbox->min[1], (mbc->c[1]) - mbc->r) && DB_GE(bbox->max[1], (mbc->c[1]) + mbc->r))
		{
			return true;
		}
	}

	return false;
}

/*does mbc covers bbox? OK*/
bool covers_from_bbox(const MBC *mbc, const BBox *bbox)
{
	double p1[2], p2[2];

	int anyBorder = 0;
	if (DB_IS_EQUAL(euclidean_distance(bbox->min, mbc->c), mbc->r))
	{
		anyBorder++;
	}
	if (DB_IS_EQUAL(euclidean_distance(bbox->max, mbc->c), mbc->r))
	{
		anyBorder++;
	}

	p1[0] = bbox->min[0];
	p1[1] = bbox->max[1];
	p2[0] = bbox->max[0];
	p2[1] = bbox->min[1];
	if (DB_IS_EQUAL(euclidean_distance(p1, mbc->c), mbc->r))
	{
		anyBorder++;
	}
	if (DB_IS_EQUAL(euclidean_distance(p2, mbc->c), mbc->r))
	{
		anyBorder++;
	}
	if (anyBorder == 0)
	{
		return false;
	}

	/* checks if the distance between bbox points and mbc center is greater than mbcs radius, if so = false*/
	if (DB_GT(euclidean_distance(bbox->min, mbc->c), mbc->r))
	{
		return false;
	}
	if (DB_GT(euclidean_distance(bbox->max, mbc->c), mbc->r))
	{
		return false;
	}
	if (DB_GT(euclidean_distance(p1, mbc->c), mbc->r))
	{
		return false;
	}
	if (DB_GT(euclidean_distance(p2, mbc->c), mbc->r))
	{
		return false;
	}
	return true;
}

/* is mbc equal to bbox? only if both are points OK*/
bool equal_from_bbox(const MBC *mbc, const BBox *bbox)
{
	/*check if mbc radius is 0 */
	if (DB_IS_EQUAL(mbc->r, 0))
	{
		/*check if bbox is a point */
		if (DB_IS_EQUAL(bbox->min[0], bbox->max[0]) && DB_IS_EQUAL(bbox->min[1], bbox->max[1]))
		{
			/* check if points match */
			if (DB_IS_EQUAL(bbox->min[0], mbc->c[0]) && DB_IS_EQUAL(bbox->min[1], mbc->c[1]))
			{
				return true;
			}
		}
	}
	return false;
}

/* is mbc inside or coveredBy bbox? OK*/
bool inside_or_coveredBy_from_bbox(const MBC *mbc, const BBox *bbox)
{
	/*checking if each side of mbc is inside bbox or on its border*/
	if (DB_LE(bbox->min[0], (mbc->c[0]) - mbc->r) && DB_GE(bbox->max[0], (mbc->c[0]) + mbc->r))
	{
		if (DB_LE(bbox->min[1], (mbc->c[1]) - mbc->r) && DB_GE(bbox->max[1], (mbc->c[1]) + mbc->r))
		{
			return true;
		}
	}

	return false;
}

/* does mbc contain or covers bbox? OK*/
bool contain_or_covers_from_bbox(const MBC *mbc, const BBox *bbox)
{
	double p1[2], p2[2];

	/* checks if the distance between bbox points and mbc center is greater than mbcs radius, if so = false*/
	if (DB_GT(euclidean_distance(bbox->min, mbc->c), mbc->r))
	{
		return false;
	}
	if (DB_GT(euclidean_distance(bbox->max, mbc->c), mbc->r))
	{
		return false;
	}

	p1[0] = bbox->min[0];
	p1[1] = bbox->max[1];
	p2[0] = bbox->max[0];
	p2[1] = bbox->min[1];
	if (DB_GT(euclidean_distance(p1, mbc->c), mbc->r))
	{
		return false;
	}
	if (DB_GT(euclidean_distance(p2, mbc->c), mbc->r))
	{
		return false;
	}
	return true;
}

/* do bbox and mbc meet? */
bool meet_from_bbox(const MBC *mbc, const BBox *bbox)
{
	int touches1 = 0, touches2 = 0;
	double p1[2], p2[2];
	/*1st check: are mbc and bbox not inside one another?*/
	if (inside_or_coveredBy_from_bbox(mbc, bbox) || contain_or_covers_from_bbox(mbc, bbox))
	{
		return false;
	}

	/*2nd check: is mbc touching how many of bbox borders?*/
	if (DB_IS_EQUAL(bbox->min[0], (mbc->c[0]) - mbc->r) && (DB_LE(bbox->min[1], mbc->c[1]) && DB_GE(bbox->max[1], mbc->c[1]))) //left: mbc center minus radius X axis / y axis in inside bbox min and max
	{
		touches1++;
	}
	if (DB_IS_EQUAL(bbox->max[0], (mbc->c[0]) + mbc->r) && (DB_LE(bbox->min[1], mbc->c[1]) && DB_GE(bbox->max[1], mbc->c[1]))) //right: mbc center plus radius X axis / y axis in inside bbox min and max
	{
		touches1++;
	}

	if (DB_IS_EQUAL(bbox->min[1], (mbc->c[1]) - mbc->r) && (DB_LE(bbox->min[0], mbc->c[0]) && DB_GE(bbox->max[0], mbc->c[0]))) //bottom: mbc center minus radius Y axis / x axis in inside bbox min and max
	{
		touches1++;
	}
	if (DB_IS_EQUAL(bbox->max[1], (mbc->c[1]) + mbc->r) && (DB_LE(bbox->min[0], mbc->c[0]) && DB_GE(bbox->max[0], mbc->c[0]))) //bottom: mbc center plus radius Y axis / x axis in inside bbox min and max
	{
		touches1++;
	}

	if (touches1 > 1)
	{
		return false;
	}

	/*3rd check: is bbox toucing how many of mbc borders?*/
	if (DB_IS_EQUAL(euclidean_distance(bbox->min, mbc->c), mbc->r))
	{
		touches2++;
	}
	if (DB_IS_EQUAL(euclidean_distance(bbox->max, mbc->c), mbc->r))
	{
		touches2++;
	}

	p1[0] = bbox->min[0];
	p1[1] = bbox->max[1];
	p2[0] = bbox->max[0];
	p2[1] = bbox->min[1];
	if (DB_IS_EQUAL(euclidean_distance(p1, mbc->c), mbc->r))
	{
		touches2++;
	}
	if (DB_IS_EQUAL(euclidean_distance(p2, mbc->c), mbc->r))
	{
		touches2++;
	}

	if (touches2 > 1)
	{
		return false;
	}

	/*if in tests 2 or 3, one touch was found = true*/
	return (touches1 == 1) || (touches2 == 1);
}

/* does mbc and bbox overlap? */
bool overlap_from_bbox(const MBC *mbc, const BBox *bbox)
{
	return !inside_or_coveredBy_from_bbox(mbc, bbox) && !contain_or_covers_from_bbox(mbc, bbox) && intersect_from_bbox(mbc, bbox) && !meet_from_bbox(mbc, bbox);
}

bool mbc_check_predicate_from_bbox(const SpatialApproximation *ap, const SpatialApproximation *bbox, uint8_t predicate)
{
	if (bbox->type == BBOX_TYPE)
	{
		MBC_APPROX *mbc_ap = (void *)ap;
		BBOX_APPROX *bbox_ap = (void *)bbox;
		MBC *mbc = &(mbc_ap->mbc);
		BBox *bbox1 = &(bbox_ap->bbox);

		switch (predicate)
		{
		case INTERSECTS:
			return intersect_from_bbox(mbc, bbox1);
		case DISJOINT:
			return !intersect_from_bbox(mbc, bbox1);
		case OVERLAP:
			return overlap_from_bbox(mbc, bbox1);
		case MEET:
			return meet_from_bbox(mbc, bbox1);
		case INSIDE:
			return inside_from_bbox(mbc, bbox1);
		case CONTAINS:
			return contains_from_bbox(mbc, bbox1);
		case COVEREDBY:
			return coveredBy_from_bbox(mbc, bbox1);
		case COVERS:
			return covers_from_bbox(mbc, bbox1);
		case EQUAL:
			return equal_from_bbox(mbc, bbox1);
		case INSIDE_OR_COVEREDBY:
			return inside_or_coveredBy_from_bbox(mbc, bbox1);
		case CONTAINS_OR_COVERS:
			return contain_or_covers_from_bbox(mbc, bbox1);
		default:
			return false;
		}
		return false;
	}
	else
	{
		_DEBUG(WARNING, "Wrong type for bbox approximation\n");
		return false;
	}
}

size_t serialize_mbc(const SpatialApproximation *ap, uint8_t *buf)
{
	uint8_t *loc;
	int i;
	MBC_APPROX *mbc_ap = (void *)ap;

	loc = buf;
	memcpy(loc, &(ap->id), sizeof(uint64_t));
	loc += sizeof(uint64_t);

	memcpy(loc, &(mbc_ap->mbc.r), sizeof(double));
	loc += sizeof(double);

	for (i = 0; i < NUM_OF_DIM; i++)
	{
		memcpy(loc, &(mbc_ap->mbc.c[i]), sizeof(double));
		loc += sizeof(double);
	}
	return (size_t)(loc - buf);
}

size_t mbc_size()
{
	return (sizeof(double) * (NUM_OF_DIM + 1)) + sizeof(uint64_t);
}

bool write_mbc_to_file(const FileSpecification *fs, const SpatialApproximation **ap, int n)
{
	size_t total_size = (n * mbc_size()) + sizeof(int), finalsize;
	uint8_t *buf, *loc;
	int i;

	if (total_size > fs->page_size)
	{
		_DEBUG(ERROR, "Total size is bigger than page size\n");
		return false;
	}

	if (fs->io_access == DIRECT_ACCESS)
	{
		if (posix_memalign((void **)&buf, fs->page_size, fs->page_size))
		{
			_DEBUG(ERROR, "Allocation failed at MBC write_to_file\n");
			return false;
		}
	}
	else
	{
		buf = lwalloc(fs->page_size);
	}

	loc = buf;

	memcpy(loc, &n, sizeof(int));
	loc += sizeof(int);

	for (i = 0; i < n; i++)
	{
		loc += serialize_mbc(ap[i], loc);
	}

	finalsize = (size_t)(loc - buf);

	if (total_size != finalsize)
	{
		_DEBUGF(WARNING, "Return size (%zu) not equal to expected size (%zu)\n", finalsize, total_size);
		return false;
	}
	append_page(fs, buf);

	if (fs->io_access == DIRECT_ACCESS)
    {
        free(buf);
    }
    else
    {
        lwfree(buf);
    }

	return true;
}

SpatialApproximation *recreate_mbc(double r, double *c, uint64_t id)
{
	static const SpatialApproximationInterface vtable = {convert_mbc_to_geom, mbc_check_predicate_from_bbox, destroy_mbc};
	static SpatialApproximation base = {&vtable};
	MBC_APPROX *approx;
	int i;

	MBC mbc;
	mbc.r = r;

	for (i = 0; i < NUM_OF_DIM; i++)
	{
		mbc.c[i] = c[i];
	}

	base.type = MBC_TYPE;
	base.id = id;
	approx = lwalloc(sizeof(MBC_APPROX));
	memcpy(&approx->base, &base, sizeof(base));
	approx->mbc = mbc;

	return &approx->base;
}

SpatialApproximation *read_serialized_mbc(uint8_t *buf, size_t *size)
{
	uint8_t *start_ptr = buf;
	double r, c[NUM_OF_DIM];
	uint64_t id;
	int i;
	SpatialApproximation *ap;

	memcpy(&id, buf, sizeof(uint64_t));
	buf += sizeof(uint64_t);

	memcpy(&r, buf, sizeof(double));
	buf += sizeof(double);

	for (i = 0; i < NUM_OF_DIM; i++)
	{
		memcpy(&(c[i]), buf, sizeof(double));
		buf += sizeof(double);
	}

	ap = recreate_mbc(r, c, id);

	*size = (size_t)(buf - start_ptr);
	return ap;
}

SpatialApproximation **read_mbc_from_file(const FileSpecification *fs, int page, int *n)
{
	uint8_t *buf, *loc;
	int i;
	size_t size;
	SpatialApproximation **ap;

	if (fs->io_access == DIRECT_ACCESS)
	{
		if (posix_memalign((void **)&buf, fs->page_size, fs->page_size))
		{
			_DEBUG(ERROR, "Allocation failed at MBC read_from_file\n");
			return NULL;
		}
	}
	else
	{
		buf = lwalloc(fs->page_size);
	}

	disk_read_one_page(fs, page, buf);

	loc = buf;

	memcpy(n, loc, sizeof(int));
	loc += sizeof(int);

	if (*n <= 0)
	{
		_DEBUGF(WARNING, "Number of approximations in page is invalid: %d\n", *n);
		return NULL;
	}

	ap = lwalloc((*n) * sizeof(SpatialApproximation *));

	for (i = 0; i < (*n); i++)
	{
		ap[i] = read_serialized_mbc(loc, &size);
		loc += size;
	}

	if (fs->io_access == DIRECT_ACCESS)
    {
        free(buf);
    }
    else
    {
        lwfree(buf);
    }

	return ap;
}

SpatialApproximation *create_mbc(LWGEOM *geom, uint64_t id)
{
	static const SpatialApproximationInterface vtable = {convert_mbc_to_geom, mbc_check_predicate_from_bbox, destroy_mbc};
	static SpatialApproximation base = {&vtable};
	MBC_APPROX *approx;

	MBC mbc;
	convert_geom_to_mbc(geom, &mbc);

	base.type = MBC_TYPE;
	base.id = id;
	approx = lwalloc(sizeof(MBC_APPROX));
	memcpy(&approx->base, &base, sizeof(base));
	approx->mbc = mbc;

	return &approx->base;
}

void destroy_mbc(SpatialApproximation *ap)
{
	MBC_APPROX *mbc_ap = (void *)ap;
	lwfree(mbc_ap);
}