#ifndef COMMAND_ANI
#define COMMAND_ANI
#include <argp.h>
#include <math.h>
#include <tgmath.h>
#include <assert.h>

#include "global_basic.h"
#include "kssdlib_sort.h"
#include "sketch_rearrange.h"
#include "model_ani.h"

typedef struct ani_opt
{
	int fmt;  // print out format: 0:detail; 1, aafD; 2. 1-ani
	double c; // minimal distance to enrolled sketches
	int p;	  // threads
	bool d;	  // diagnal
	bool v;	  // naive model ?
	int e;
	float afcut;
	float anicut;
	char index[PATHLEN];
	char qrydir[PATHLEN];
	char refdir[PATHLEN];
	char outf[PATHLEN];
	char gl[PATHLEN]; // genome list with selection code
	char model[PATHLEN];
	int num_remaining_args;
	char **remaining_args;
} ani_opt_t;

typedef struct
{
	uint32_t num_ctx; // including with confilict obj
	uint32_t num_conflictobj;
} num_ctx_cfltobj_t;

typedef struct
{
	uint64_t arrlen;
	uint32_t num_ctx;
	uint32_t num_conflictobj;
	double numgids_perctx;
	int infile_num;
	num_ctx_cfltobj_t *num_ctx_cfltobj_arr;
} sort_sketch_summary_t;

typedef struct
{
	uint32_t id;
	float ani;
} idani_t;

typedef struct
{
	uint32_t diff_obj;
	uint32_t diff_obj_section;
} obj_section_t;

typedef struct
{
	uint32_t num_ctx;
	uint32_t num_mut2_ctx;
} ctx_mut2_t;

int cmd_ani(struct argp_state *state);
int compute_ani(ani_opt_t *ani_opt);
int mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt);
int sparse_mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt);
int compare_u32(const void *a, const void *b);
void sort_arrays(uint64_t *a, uint32_t *b, uint32_t n);
ctxgidobj_t *comb_sortedsketch64_2sortedcomb_ctxgid64obj32(unify_sketch_t *result);
sort_sketch_summary_t *summarize_ctxgidobj_arr(ctxgidobj_t *ctxgidobj_arr, uint64_t *sketch_index, uint32_t arrlen, int infile_num);
void free_sort_sketch_summary(sort_sketch_summary_t *sort_sketch_summary);
// void sorted_ctxgidobj_arr2triangle (ctxgidobj_t* ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary);
void sorted_ctxgidobj_arrXcomb_sortedsketch64(unify_sketch_t *qry_result, ctxgidobj_t *ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary);
// void comb_sortedsketch64Xcomb_sortedsketch64 ( unify_sketch_t* ref_result, unify_sketch_t* qry_result );
void comb_sortedsketch64Xcomb_sortedsketch64(ani_opt_t *ani_opt);

size_t *find_first_occurrences_AT_ctxgidobj_arr(const uint64_t *a, size_t a_size, const ctxgidobj_t *b, size_t b_size);
// void get_ani_features_from_two_sorted_ctxobj64 (const uint64_t *a, size_t n,  const uint64_t *b, size_t m, ani_features_t* ani_features);
// print functions
void ani_block_print_matrix(int ref_infile_num, int qry_gid_offset, int this_block_size, uint64_t *qry_sketch_index, ctx_mut2_t *ctx, obj_section_t *obj, char (*qryfname)[PATHLEN], FILE *outfp, ani_opt_t *ani_opt);
void ani_block_print(int ref_infile_num, int qry_gid_offset, int this_block_size, uint64_t *ref_sketch_index, uint64_t *qry_sketch_index, ctx_mut2_t *ctx, obj_section_t *obj, char (*refname)[PATHLEN], char (*qryfname)[PATHLEN], uint32_t *num_passid_block, idani_t **sort_idani_block, FILE *outfp);

// inline functions
// 1. 3-way linear model distance (see model_ani.h): static inline double lm3ways_dist_from_features(ani_features_t *features)
// 2. naive distance: for where 3-way linear model is not applicable e.g. unassembled genomes , or Eukaryotic genomes?
static inline double get_naive_dist(ani_features_t *features)
{
	if (features->XnY_ctx == 0)
		return 1;
	double ratio = (double)(features->N_diff_obj_section + EPSILON) / (features->N_diff_obj + EPSILON);
	double dist0 = (double)features->N_diff_obj / (features->XnY_ctx + features->N_diff_obj);
	double final_dist = 1 - pow((1 - dist0), ratio);
	return final_dist / 7; // 7 is experically determined;
}
// 3. generic distance function from 1. or 2.
typedef double (*get_generic_dist_from_features_fn)(ani_features_t *features);
extern get_generic_dist_from_features_fn get_generic_dist_from_features ;
// called by comb_sortedsketch64Xcomb_sortedsketch64()
// cautions: a and b must be obj conflition removed arrays.
static inline void get_ani_features_from_two_sorted_ctxobj64(const uint64_t *a, size_t n,
															 const uint64_t *b, size_t m, ani_features_t *ani_features)
{
	uint8_t nobjbits = Bitslen.obj;
	uint64_t objmask = (1UL << nobjbits) - 1;
	size_t i = 0, j = 0;
	memset(ani_features, 0, sizeof(ani_features_t));

	while (i < n && j < m)
	{

		if (a[i] >> nobjbits == b[j] >> nobjbits)
		{
			ani_features->XnY_ctx++;
			uint32_t has_diff_obj = (uint32_t)(a[i] & objmask) ^ (b[j] & objmask);
			if (has_diff_obj)
			{
				ani_features->N_diff_obj++;
				int num_diff_obj_section = 0;
				for (int k = 0; k < nobjbits / 2; k++)
				{
					if (has_diff_obj & (3U << (2 * k)))
						num_diff_obj_section++;
				}
				ani_features->N_diff_obj_section += num_diff_obj_section;
				if (num_diff_obj_section > 1)
					ani_features->N_mut2_ctx++;
			}
			i++;
			j++;
		}
		else if (a[i] < b[j])
		{ // else if (a[i] >> nobjbits < b[j] >> nobjbits)
			i++;
		}
		else
		{
			j++;
		}
	}
}

static int compare_idani_desc(const void *a, const void *b)
{
	const idani_t *itemA = (const idani_t *)a;
	const idani_t *itemB = (const idani_t *)b;
	return (itemA->ani < itemB->ani) - (itemA->ani > itemB->ani);
}

#define MCTX(L, X, Y) (ctx[(size_t)((L) * (X) + (Y))])
#define MOBJ(L, X, Y) (obj[(size_t)((L) * (X) + (Y))])
static inline void count_ctx_obj_frm_comb_sketch_section(ctx_mut2_t *ctx, obj_section_t *obj, ctxgidobj_t *ctxgidobj_arr, size_t ref_sksize, int ref_gnum, int section_gnum, uint64_t *section_sk, uint64_t *section_skidx, uint32_t *num_passid_block, idani_t **sort_idani_block, ani_opt_t *ani_opt)
{
	//	printf("%d\n",Bitslen.obj/2); exit(0);
	uint64_t gidmask = UINT64_MAX >> (64 - GID_NBITS), objmask = (1UL << Bitslen.obj) - 1;
#pragma omp parallel for num_threads(ani_opt->p) schedule(guided)
	for (int i = 0; i < section_gnum; i++)
	{
		uint64_t *a = section_sk + (section_skidx[i] - section_skidx[0]);
		size_t a_size = section_skidx[i + 1] - section_skidx[i];
		assert(a_size > 0);
		size_t *idx = find_first_occurrences_AT_ctxgidobj_arr(a, a_size, ctxgidobj_arr, ref_sksize);

		for (int j = 0; j < a_size; j++)
		{
			if (idx[j] == SIZE_MAX)
				continue;
			// Skip when no findings, or conflict objects (adjacent elements with the same context)
			// if ((j > 0 && (a[j] >> Bitslen.obj) == (a[j - 1] >> Bitslen.obj)) ||
			//	(j < a_size - 1 && (a[j] >> Bitslen.obj) == (a[j + 1] >> Bitslen.obj)))
			//	continue;
			for (int d = idx[j];; d++)
			{
				if ((ctxgidobj_arr[d].ctxgid >> Bitslen.gid) != (a[j] >> Bitslen.obj))
					break;
				uint32_t gid = ctxgidobj_arr[d].ctxgid & gidmask;
				MCTX(ref_gnum, i, gid).num_ctx++;
				uint32_t has_diff_obj = (uint32_t)(a[j] & objmask) ^ ctxgidobj_arr[d].obj;
				if (has_diff_obj)
				{
					MOBJ(ref_gnum, i, gid).diff_obj++;
					// count diff_obj_section
					int num_diff_obj_section = 0;
					for (int k = 0; k < Bitslen.obj / 2; k++)
					{
						if (has_diff_obj & (3U << (2 * k)))
							num_diff_obj_section++;
					}
					MOBJ(ref_gnum, i, gid).diff_obj_section += num_diff_obj_section;
					if (num_diff_obj_section > 1)
						MCTX(ref_gnum, i, gid).num_mut2_ctx++;
				}
			}
		}
		free(idx);
		if (ani_opt->fmt == 0)
		{ // print details  >0:matrix no need to sort
			// sorting ref gid by ani descendingly
			num_passid_block[i] = 0;
			for (int j = 0; j < ref_gnum; j++)
			{
				if ((float)MCTX(ref_gnum, i, j).num_ctx / a_size < ani_opt->afcut)
					continue;
				float dist = (float)MOBJ(ref_gnum, i, j).diff_obj / MCTX(ref_gnum, i, j).num_ctx;
				float ani = 1 - dist;
				if (ani < ani_opt->anicut)
					continue;
				sort_idani_block[i][num_passid_block[i]].id = j;
				sort_idani_block[i][num_passid_block[i]].ani = ani;
				num_passid_block[i]++;
			}
			qsort(sort_idani_block[i], num_passid_block[i], sizeof(idani_t), compare_idani_desc);
		}
	}
}

#endif
