#include <string.h>
#include <math.h>
#include "SpatialApproximation_utils.h"
#include "../main/log_messages.h"

SpatialApproximation **spatialapproximation_builder(int total, LWGEOM **geoms, uint8_t type, int *ids)
{
	SpatialApproximation **approx;
	int i;
	SpatialApproximation *(*create)(LWGEOM *, uint64_t);
	if (total < 1 || geoms == NULL)
		return NULL;

	approx = lwalloc(total * sizeof(SpatialApproximation *));

	create = &create_bbox_approx;
	if (type == MBC_TYPE)
	{
		create = &create_mbc;
	}
	else if (type == RMBR_TYPE)
	{
		create = &create_rmbr;
	}
	else if (type == RMBP_TYPE)
	{
		create = &create_rmbp;
	}
	else if (type == MBE_TYPE)
	{
		create = &create_mbe;
	}
	else if ((type & 0b1111) == N_CORNER_TYPE)
	{
		create = NULL;
	}

	for (i = 0; i < total; ++i)
	{
		if (create)
		{
			approx[i] = create(geoms[i], ids[i]);
		}
		else
		{
			approx[i] = create_n_corner(geoms[i], type, ids[i]);
		}
	}

	return approx;
}

SpatialApproximation **spatialapproximation_get_from_list(const int *ids, int n, uint8_t type, const FileSpecification *fs)
{
	SpatialApproximation **ap;
	size_t ap_size, size;
	int num_per_page, i, page_no, no;
	uint8_t *buf, *loc;
	if (n <= 0)
	{
		_DEBUG(ERROR, "Invalid list size\n");
		return NULL;
	}
	ap = lwalloc(n * sizeof(SpatialApproximation *));
	ap_size = spatialapproximation_size(type);
	num_per_page = floor(((double)fs->page_size - sizeof(uint32_t)) / (double)ap_size);

	for (i = 0; i < n; ++i)
	{
		//Calculating page index
		page_no = floor((double)ids[i] / (double)num_per_page);

		if ((ids[i] % num_per_page) == 0)
		{
			page_no--;
		}

		//Reading page
		if (fs->io_access == DIRECT_ACCESS)
		{
			if (posix_memalign((void **)&buf, fs->page_size, fs->page_size))
			{
				_DEBUG(ERROR, "Allocation failed at spatialapproximation_get_from_list\n");
				return NULL;
			}
		}
		else
		{
			buf = lwalloc(fs->page_size);
		}

		disk_read_one_page(fs, page_no, buf);

		loc = buf;

		memcpy(&no, loc, sizeof(int));
		loc += sizeof(int);

		if (no <= 0)
		{
			_DEBUGF(WARNING, "Number of approximations in page is invalid: %d\n", no);
			return NULL;
		}

		//Positioning buffer and getting approximation
		loc += (ids[i] - (page_no * num_per_page) - 1) * ap_size;

		ap[i] = spatialapproximation_read_from_serialization(type, loc, &size);

		//Memory management
		if (fs->io_access == DIRECT_ACCESS)
		{
			free(buf);
		}
		else
		{
			lwfree(buf);
		}
	}

	return ap;
}

int *spatialapproximation_filter(const SpatialApproximation *bbox, int n, SpatialApproximation **ap, uint8_t predicate, int *m)
{
	if (bbox->type == BBOX_TYPE)
	{
		int *temp = lwalloc(n * sizeof(uint64_t)), *ret;
		int c = 0, i;

		for (i = 0; i < n; ++i)
		{
			if (spatialapproximation_check_predicate(ap[i], bbox, predicate))
			{
				temp[c] = ap[i]->id;
				c++;
			}
		}

		*m = c;

		ret = lwalloc(c * sizeof(uint64_t));

		for (i = 0; i < c; i++)
		{
			ret[i] = temp[i];
		}

		lwfree(temp);
		return ret;
	}
	else
	{
		_DEBUG(ERROR, "BBox Spatial Approximation received is not BBOX_TYPE\n");
		return NULL;
	}
}

void spatialapproximation_free_array(int n, SpatialApproximation **ap)
{
	int i;
	for (i = 0; i < n; i++)
	{
		spatialapproximation_free(ap[i]);
	}
	lwfree(ap);
}

char *get_approx_name(uint8_t type)
{
	if (type == MBC_TYPE)
	{
		return "mbc";
	}
	else if (type == RMBR_TYPE)
	{
		return "rmbr";
	}
	else if (type == RMBP_TYPE)
	{
		return "rmbp";
	}
	else if (type == MBE_TYPE)
	{
		return "mbe";
	}
	else if (type == N_CORNER_5_TYPE)
	{
		return "5corner";
	}
	else if (type == N_CORNER_4_TYPE)
	{
		return "4corner";
	}
	else
	{
		return "bbox";
	}
}

uint8_t get_approx_type(char *file_path)
{
	char *last_dot = strrchr(file_path, '.');

	if (strcmp(last_dot + 1, "mbc") == 0)
	{
		return MBC_TYPE;
	}
	else if (strcmp(last_dot + 1, "rmbr") == 0)
	{
		return RMBR_TYPE;
	}
	else if (strcmp(last_dot + 1, "rmbp") == 0)
	{
		return RMBP_TYPE;
	}
	else if (strcmp(last_dot + 1, "mbe") == 0)
	{
		return MBE_TYPE;
	}
	else if (strcmp(last_dot + 1, "4corner") == 0)
	{
		return N_CORNER_4_TYPE;
	}
	else if (strcmp(last_dot + 1, "5corner") == 0)
	{
		return N_CORNER_5_TYPE;
	}
	else
	{
		return BBOX_TYPE;
	}
}