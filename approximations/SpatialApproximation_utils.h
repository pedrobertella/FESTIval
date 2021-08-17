#ifndef _SPATIAL_APPROXIMATION_UTILS_H
#define _SPATIAL_APPROXIMATION_UTILS_H

#include "SpatialApproximation.h"
#include "../main/io_handler.h"
#include "rmbp_handler.h"
#include "rmbr_handler.h"
#include "n-corner_handler.h"
#include "mbc_handler.h"
#include "mbe_handler.h"
#include "bbox_approx_handler.h"

/* Utils */

static inline size_t spatialapproximation_size(const uint8_t type)
{
	if (type == BBOX_TYPE)
	{
		return bbox_size();
	}
	else if (type == MBC_TYPE)
	{
		return mbc_size();
	}
	else if (type == RMBR_TYPE)
	{
		return rmbr_size();
	}
	else if (type == RMBP_TYPE)
	{
		return rmbp_size();
	}
	else if (type == MBE_TYPE)
	{
		return mbe_size();
	}
	else if (type == N_CORNER_4_TYPE)
	{
		return n_4_corner_size();
	}
	else if (type == N_CORNER_5_TYPE)
	{
		return n_5_corner_size();
	}
	return 0;
}

extern SpatialApproximation **spatialapproximation_builder(int total, LWGEOM **geoms, uint8_t type, int *ids);

extern SpatialApproximation **spatialapproximation_get_from_list(const int *ids, int n, uint8_t type, const FileSpecification *fs);

extern int *spatialapproximation_filter(const SpatialApproximation *bbox, int n, SpatialApproximation **ap, uint8_t predicate, int *m);

extern void spatialapproximation_free_array(int n, SpatialApproximation **ap);

extern char *get_approx_name(uint8_t type);

extern uint8_t get_approx_type(char *file_path);

/* I/O functions */

static inline bool spatialapproximation_write_to_file(uint8_t type, const FileSpecification *fs, const SpatialApproximation **ap, int n)
{
	switch (type)
	{
	case RMBP_TYPE:
		return write_rmbp_to_file(fs, ap, n);
	case RMBR_TYPE:
		return write_rmbr_to_file(fs, ap, n);
	case MBC_TYPE:
		return write_mbc_to_file(fs, ap, n);
	case BBOX_TYPE:
		return write_bbox_to_file(fs, ap, n);
	case MBE_TYPE:
		return write_mbe_to_file(fs, ap, n);
	case N_CORNER_4_TYPE:
	case N_CORNER_5_TYPE:	
		return write_n_corner_to_file(fs, ap, n);
	default:
		break;
	}
	return false;
}

static inline SpatialApproximation **spatialapproximation_read_from_file(uint8_t type, const FileSpecification *fs, int page, int *n)
{
	switch (type)
	{
	case RMBP_TYPE:
		return read_rmbp_from_file(fs, page, n);
	case RMBR_TYPE:
		return read_rmbr_from_file(fs, page, n);
	case MBC_TYPE:
		return read_mbc_from_file(fs, page, n);
	case BBOX_TYPE:
		return read_bbox_from_file(fs, page, n);
	case MBE_TYPE:
		return read_mbe_from_file(fs, page, n);
	case N_CORNER_4_TYPE:
	case N_CORNER_5_TYPE:	
		return read_n_corner_from_file(fs, page, n);
	default:
		break;
	}
	return NULL;
}

static inline SpatialApproximation *spatialapproximation_read_from_serialization(uint8_t type, uint8_t *buf, size_t *size)
{
	switch (type)
	{
	case RMBP_TYPE:
		return read_serialized_rmbp(buf, size);
	case RMBR_TYPE:
		return read_serialized_rmbr(buf, size);
	case MBC_TYPE:
		return read_serialized_mbc(buf, size);
	case BBOX_TYPE:
		return read_serialized_bbox(buf, size);
	case MBE_TYPE:
		return read_serialized_mbe(buf, size);
	case N_CORNER_4_TYPE:
	case N_CORNER_5_TYPE:	
		return read_serialized_n_corner(buf, size);
	default:
		break;
	}
	return NULL;
}

#endif /* _SPATIAL_APPROXIMATION_UTILS_H */