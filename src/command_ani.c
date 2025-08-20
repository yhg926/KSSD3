#include "command_ani.h"
#include "command_matrix.h"
#include "global_basic.h"
#include "kssdlib_sort.h"
#include "sketch_rearrange.h"
// #include "command_sketch.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <err.h>
#include <errno.h>
#include <math.h>
#include <libgen.h>
#include <dirent.h>
#include <omp.h>
#include <stdatomic.h>
#include <ctype.h>
#include "../klib/kstring.h" // from klib
#include "../klib/khash.h"
// #include "../klib/khashl.h"
#define GID_NBITS 20 // 2^20, 1M
#define CONFLICT_OBJ UINT32_MAX
// pulic vars
const char gid_obj_prefix[] = "gidobj", ctx_idx_prefix[] = "ctx.index";
extern const char sorted_comb_ctxgid64obj32[];
extern double C9O7_98[6], C9O7_96[6];
size_t file_size;

const char print_header[] = "Qry\tRef\tXnY_ctx\tQry_align_fraction\tRef_align_fraction\tN_diff_obj\tN_diff_obj_section\tN_mut2_ctx\tANI";
const char select_metrics_header[5][15] = {"StDist", "MashD", "AafD", "MashD_if_far", "AafD_if_far"};
double select_metrics_dist[5] = {1, 1, 1, 1, 1}; // default values

#define BLOCK_SIZE (4096) // #of qry genomes per batch, for mem_eff handling
int mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt)
{

	dim_sketch_stat_t *ref_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->refdir, sketch_stat), &file_size);
	int ref_infile_num = ref_dim_sketch_stat->infile_num;
	// read index
	size_t ctxgidobj_arr_fsize;
	uint64_t *ref_sketch_index = read_from_file(test_get_fullpath(ani_opt->refdir, idx_sketch_suffix), &file_size);
	assert(file_size == (ref_infile_num + 1) * sizeof(ref_sketch_index[0]));
	size_t ref_sketch_size = ref_sketch_index[ref_infile_num];
	ctxgidobj_t *sortedcomb_ctxgid64obj32 = read_from_file(test_get_fullpath(ani_opt->refdir, sorted_comb_ctxgid64obj32), &ctxgidobj_arr_fsize);
	assert(ctxgidobj_arr_fsize == ref_sketch_size * sizeof(sortedcomb_ctxgid64obj32[0]));

	dim_sketch_stat_t *qry_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->qrydir, sketch_stat), &file_size);
	int qry_infile_num = qry_dim_sketch_stat->infile_num;
	assert(qry_dim_sketch_stat->hash_id == ref_dim_sketch_stat->hash_id);
	uint64_t *qry_sketch_index = read_from_file(test_get_fullpath(ani_opt->qrydir, idx_sketch_suffix), &file_size);
	size_t qry_sketch_size = qry_sketch_index[qry_infile_num];

	int block_size = BLOCK_SIZE;
	int offset_gid = 0;
	FILE *fp;
	assert((fp = fopen(test_get_fullpath(ani_opt->qrydir, combined_sketch_suffix), "rb")) != NULL);
	uint64_t *tmp_ctxobj = malloc(qry_sketch_size * sizeof(uint64_t));
	// uint32_t *ctx = malloc(ref_infile_num * block_size * sizeof(uint32_t));
	// uint32_t *obj = malloc( ref_infile_num * block_size * sizeof(uint32_t))  ;
	ctx_mut2_t *ctx = malloc(ref_infile_num * block_size * sizeof(ctx_mut2_t));
	obj_section_t *obj = malloc(ref_infile_num * block_size * sizeof(obj_section_t));
	// for order id by descending ani
	uint32_t *num_passid_block = malloc(block_size * sizeof(uint32_t));
	idani_t **sort_idani_block = malloc(block_size * sizeof(idani_t *));
	for (int i = 0; i < block_size; i++)
		sort_idani_block[i] = malloc(ref_infile_num * sizeof(idani_t));

	char (*refname)[PATHLEN] = (char (*)[PATHLEN])(ref_dim_sketch_stat + 1);
	char (*qryname)[PATHLEN] = (char (*)[PATHLEN])(qry_dim_sketch_stat + 1);

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");

	/* load model
	if (ani_opt->model[0] == '\0')
		err(EXIT_FAILURE, "%s(): need specify model file using -M ", __func__);
	init_model(ani_opt->model); // f8C9O7_model.xgb
	*/

	// printf header
	if (ani_opt->fmt)
	{ // matrix format
		for (int i = 0; i < ref_infile_num; i++)
			fprintf(outfp, "\t%s", refname[i]);
		fprintf(outfp, "\n");
	}
	else
		fprintf(outfp, "%s\t%s\n", print_header, select_metrics_header[abs(ani_opt->s) - 1]);

	for (int b = 0; b <= qry_infile_num / block_size; b++)
	{

		int this_block_size = (b == qry_infile_num / block_size) ? (qry_infile_num % block_size) : block_size;
		int this_sketch_size = qry_sketch_index[offset_gid + this_block_size] - qry_sketch_index[offset_gid];
		int read_sketch_size = fread(tmp_ctxobj, sizeof(uint64_t), this_sketch_size, fp);
		uint64_t *this_sketch_index = qry_sketch_index + offset_gid;
		assert(this_sketch_size == read_sketch_size);

		memset(ctx, 0, ref_infile_num * block_size * sizeof(ctx_mut2_t));
		memset(obj, 0, ref_infile_num * block_size * sizeof(obj_section_t)); // memset(obj,0,ref_infile_num * block_size * sizeof(uint32_t));
		count_ctx_obj_frm_comb_sketch_section(ctx, obj, sortedcomb_ctxgid64obj32, ref_sketch_size, ref_infile_num, this_block_size, tmp_ctxobj, this_sketch_index, num_passid_block, sort_idani_block, ani_opt);
		ani_block_print(ref_infile_num, b * block_size, this_block_size, ref_sketch_index, qry_sketch_index, ctx, obj, refname, qryname, num_passid_block, sort_idani_block, outfp, ani_opt, ani_opt->fmt);

		offset_gid += this_block_size;
	}

	for (int i = 0; i < block_size; i++)
		free(sort_idani_block[i]);
	free_all(ref_dim_sketch_stat, ref_sketch_index, qry_dim_sketch_stat, qry_sketch_index, tmp_ctxobj, ctx, obj, num_passid_block, sort_idani_block, NULL);
	free_read_from_file(sortedcomb_ctxgid64obj32, ctxgidobj_arr_fsize);
	fclose(fp);
	fclose(outfp);
	// clean xgb model
	// cleanup_model();
	return ctxgidobj_arr_fsize;
}

//
#define DIFF_OBJ_BITS 1
size_t dedup_with_ctxobj_counts(uint32_t *arr, size_t n, co_distance_t **ctxobj_cnt)
{
	*ctxobj_cnt = NULL;
	if (n == 0)
		return 0;

	co_distance_t *tmp_ctxobj_cnt = malloc(n * sizeof(co_distance_t));
	if (!tmp_ctxobj_cnt)
		err(EXIT_FAILURE, "%s(): tmp_ctxobj_cnt malloc failure", __func__);

	size_t j = 0;
	tmp_ctxobj_cnt[0].ctx_ct = 1;
	tmp_ctxobj_cnt[0].diff_obj = arr[0] % 2;
	arr[0] >>= DIFF_OBJ_BITS;

	for (size_t i = 1; i < n; i++)
	{
		if ((arr[i] >> DIFF_OBJ_BITS) == arr[j])
		{
			tmp_ctxobj_cnt[j].ctx_ct++;
			tmp_ctxobj_cnt[j].diff_obj += (arr[i] % 2);
		}
		else
		{
			j++;
			tmp_ctxobj_cnt[j].ctx_ct = 1;
			tmp_ctxobj_cnt[j].diff_obj = arr[i] % 2;
			arr[j] = arr[i] >> DIFF_OBJ_BITS;
		}
	}
	// Trim arrays
	co_distance_t *ctxobj_tmp = realloc(tmp_ctxobj_cnt, (j + 1) * sizeof(co_distance_t));
	*ctxobj_cnt = ctxobj_tmp ? ctxobj_tmp : tmp_ctxobj_cnt;

	return j + 1;
}

ctxgidobj_t *comb_sortedsketch64_2sortedcomb_ctxgid64obj32(unify_sketch_t *ref_result)
{
	// const_comask_init(&ref_result->stats.lco_stat_val);
	uint64_t sketch_size = ref_result->sketch_index[ref_result->infile_num];
	if (sketch_size > (float)UINT32_MAX * LD_FCTR)
		err(EXIT_FAILURE, "%s():sketch_index maximun %lu exceed UINT32_MAX*LF;%f", __func__, sketch_size, (float)UINT32_MAX * LD_FCTR);
	if (ref_result->infile_num >= (1 << GID_NBITS))
		err(EXIT_FAILURE, "%s(): genome numer %d exceed maximum:%u", __func__, ref_result->infile_num, 1 << GID_NBITS);
	if (GID_NBITS + 4 * hclen > 64)
		err(EXIT_FAILURE, "%s(): context_bits_len(%d)+gid_bits_len(%d) exceed 64", __func__, 4 * hclen, GID_NBITS);
	ctxgidobj_t *ctxgidobj_arr = ctxobj64_2ctxgidobj(ref_result->sketch_index, ref_result->comb_sketch, ref_result->infile_num, sketch_size);
	ctxgidobj_sort_array(ctxgidobj_arr, sketch_size);
	return ctxgidobj_arr;
}

sort_sketch_summary_t *summarize_ctxgidobj_arr(ctxgidobj_t *ctxgidobj_arr, uint64_t *sketch_index, uint32_t arrlen, int infile_num)
{

	uint64_t gidmask = UINT64_MAX >> (64 - GID_NBITS);
	num_ctx_cfltobj_t *num_ctx_cfltobj_arr = malloc(infile_num * sizeof(num_ctx_cfltobj_t));
	for (int i = 0; i < infile_num; i++)
	{
		num_ctx_cfltobj_arr[i].num_ctx = sketch_index[i + 1] - sketch_index[i];
		num_ctx_cfltobj_arr[i].num_conflictobj = 0;
	}

	uint32_t num_ctx = 1;
	for (uint32_t i = 1; i < arrlen; i++)
	{
		if (ctxgidobj_arr[i].ctxgid >> GID_NBITS != ctxgidobj_arr[i - 1].ctxgid >> GID_NBITS)
			num_ctx++;
		if (ctxgidobj_arr[i].ctxgid == ctxgidobj_arr[i - 1].ctxgid)
		{
			uint32_t gid = (uint32_t)ctxgidobj_arr[i].ctxgid & gidmask;
			if (num_ctx_cfltobj_arr[gid].num_ctx == (sketch_index[gid + 1] - sketch_index[gid]))
				num_ctx_cfltobj_arr[gid].num_conflictobj++;

			num_ctx_cfltobj_arr[gid].num_ctx--;
		}
	}
	//	for(int i = 0 ; i < infile_num; i++) num_conflictobj +=  num_ctx_cfltobj_arr[i].num_conflictobj;
	sort_sketch_summary_t *sort_sketch_summary = malloc(sizeof(sort_sketch_summary_t));
	sort_sketch_summary->arrlen = arrlen;
	sort_sketch_summary->numgids_perctx = (double)arrlen / num_ctx;
	sort_sketch_summary->num_ctx = num_ctx;
	sort_sketch_summary->infile_num = infile_num;
	sort_sketch_summary->num_ctx_cfltobj_arr = num_ctx_cfltobj_arr;
	return sort_sketch_summary;
}

void free_sort_sketch_summary(sort_sketch_summary_t *sort_sketch_summary)
{
	free(sort_sketch_summary->num_ctx_cfltobj_arr);
	free(sort_sketch_summary);
}

#define CTX(X, Y) (ctx[(size_t)(((X) * ((X) + 1)) / 2 + (Y))])
#define OBJ(X, Y) (obj[(size_t)(((X) * ((X) + 1)) / 2 + (Y))])

/* use sketch variants to caculate distance */
// 1.global sorted comb_sketch64 (i.e. sorted_ctxgidobj_arr or inverted index):
// slow when dist matrix is large, low memory cache efficient, but may be very fast when matrix is small?
void sorted_ctxgidobj_arr2triangle(ctxgidobj_t *ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary)
{
	uint64_t gidmask = UINT64_MAX >> (64 - GID_NBITS);
	uint32_t arrlen = sort_sketch_summary->arrlen;
	int infile_num = sort_sketch_summary->infile_num;
	uint16_t *ctx = calloc((size_t)infile_num * (infile_num + 1) / 2, sizeof(uint16_t));
	uint16_t *obj = calloc((size_t)infile_num * (infile_num + 1) / 2, sizeof(uint16_t));

	for (uint32_t i = 0, j; i < arrlen - 1; i = j)
	{
		j = i + 1;
		// find range i..j;
		for (; j < arrlen && (ctxgidobj_arr[i].ctxgid >> GID_NBITS == ctxgidobj_arr[j].ctxgid >> GID_NBITS); j++)
			;

#pragma omp parallel for num_threads(32) schedule(guided)
		for (uint32_t a = i + 1; a < j; a++)
		{
			if (ctxgidobj_arr[a].ctxgid == ctxgidobj_arr[a - 1].ctxgid || ctxgidobj_arr[a].ctxgid == ctxgidobj_arr[a + 1].ctxgid)
				continue; // with confclit object
			uint32_t x = ctxgidobj_arr[a].ctxgid & gidmask;
#pragma omp parallel for num_threads(32) schedule(guided)
			for (uint32_t b = i; b < a; b++)
			{

				if ((b > 0 && ctxgidobj_arr[b].ctxgid == ctxgidobj_arr[b - 1].ctxgid) || ctxgidobj_arr[b].ctxgid == ctxgidobj_arr[b + 1].ctxgid)
					continue;
				uint32_t y = ctxgidobj_arr[b].ctxgid & gidmask;
				CTX(x, y)
				++;
				if (ctxgidobj_arr[a].obj != ctxgidobj_arr[b].obj)
					OBJ(x, y)
				++;
			}
		}
		printf("\ri=%d", i);
	}
	for (int x = 1; x < infile_num; x++)
	{

		for (int y = 0; y < x; y++)
		{
			if (CTX(x, y) > 0)
			{
				printf("%d\t%d\t%d\t%d\t%f\n", x, y, CTX(x, y), OBJ(x, y), (float)OBJ(x, y) / CTX(x, y));
			}
		}
	}
}

// sparse_mem_eff.. seems slower than mem_eff,
int sparse_mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt)
{
	// initialize
	dim_sketch_stat_t *ref_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->refdir, sketch_stat), &file_size);
	const_comask_init(ref_dim_sketch_stat);

	uint64_t gidmask = UINT64_MAX >> (64 - GID_NBITS);
	uint64_t objmask = (1UL << Bitslen.obj) - 1;

	int ref_infile_num = ref_dim_sketch_stat->infile_num;
	// read index
	size_t ctxgidobj_arr_fsize;
	uint64_t *ref_sketch_index = read_from_file(test_get_fullpath(ani_opt->refdir, idx_sketch_suffix), &file_size);
	assert(file_size == (ref_infile_num + 1) * sizeof(ref_sketch_index[0]));
	size_t ref_sketch_size = ref_sketch_index[ref_infile_num];
	ctxgidobj_t *sortedcomb_ctxgid64obj32 = read_from_file(test_get_fullpath(ani_opt->refdir, sorted_comb_ctxgid64obj32), &ctxgidobj_arr_fsize);
	assert(ctxgidobj_arr_fsize == ref_sketch_size * sizeof(sortedcomb_ctxgid64obj32[0]));

	dim_sketch_stat_t *qry_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->qrydir, sketch_stat), &file_size);
	int qry_infile_num = qry_dim_sketch_stat->infile_num;
	assert(qry_dim_sketch_stat->hash_id == ref_dim_sketch_stat->hash_id);
	uint64_t *qry_sketch_index = read_from_file(test_get_fullpath(ani_opt->qrydir, idx_sketch_suffix), &file_size);
	size_t qry_sketch_size = qry_sketch_index[qry_infile_num];

	int block_size = BLOCK_SIZE;
	int offset_gid = 0;
	FILE *fp;
	assert((fp = fopen(test_get_fullpath(ani_opt->qrydir, combined_sketch_suffix), "rb")) != NULL);
	uint64_t *tmp_ctxobj = malloc(qry_sketch_size * sizeof(uint64_t));
	// arry of dynamic gids array per query, note gid <<= DIFF_OBJ_BITS + (ref_obj==qry_obj?0:1);  >>sparse only
	//    Vector *ref_gids_perqry_arr  = malloc( block_size * sizeof(Vector));
	// ref_gids_perqry_arr[i][0] is capacity, ref_gids_perqry_arr[i][1] is length of the inner array
	uint32_t **ref_gids_perqry_arr = malloc(block_size * sizeof(uint32_t *));
	co_distance_t **ctxobj_cnt_perqry_arr = malloc(block_size * sizeof(co_distance_t *));
	size_t *lens = malloc(block_size * sizeof(size_t));
	for (int i = 0; i < block_size; i++)
	{
		ref_gids_perqry_arr[i] = malloc(6 * sizeof(uint32_t)); // 4 for innitialized capacity
		ref_gids_perqry_arr[i][0] = 4;
		ref_gids_perqry_arr[i][1] = 0;
	}

	for (int b = 0; b <= qry_infile_num / block_size; b++)
	{

		int this_block_size = (b == qry_infile_num / block_size) ? (qry_infile_num % block_size) : block_size;
		int this_sketch_size = qry_sketch_index[offset_gid + this_block_size] - qry_sketch_index[offset_gid];
		int read_sketch_size = fread(tmp_ctxobj, sizeof(uint64_t), this_sketch_size, fp);
		uint64_t *this_sketch_index = qry_sketch_index + offset_gid;
		assert(this_sketch_size == read_sketch_size);

#pragma omp parallel for num_threads(ani_opt->p) schedule(guided)
		for (int i = 0; i < this_block_size; i++)
		{
			uint64_t *a = tmp_ctxobj + (this_sketch_index[i] - this_sketch_index[0]);
			size_t a_size = this_sketch_index[i + 1] - this_sketch_index[i];
			size_t *idx = find_first_occurrences_AT_ctxgidobj_arr(a, a_size, sortedcomb_ctxgid64obj32, ref_sketch_size);

			for (int j = 0; j < a_size; j++)
			{
				if (idx[j] == SIZE_MAX)
					continue;
				// skip conlict object;
				if ((j > 0) && ((a[j] >> Bitslen.obj) == (a[j - 1] >> Bitslen.obj)))
					continue;
				if ((j < a_size - 1) && ((a[j] >> Bitslen.obj) == (a[j + 1] >> Bitslen.obj)))
					continue;

				for (int d = idx[j];; d++)
				{
					if ((sortedcomb_ctxgid64obj32[d].ctxgid >> Bitslen.gid) != (a[j] >> Bitslen.obj))
						break;
					uint32_t gid01 = sortedcomb_ctxgid64obj32[d].ctxgid & gidmask << DIFF_OBJ_BITS;
					if ((a[j] & objmask) != sortedcomb_ctxgid64obj32[d].obj)
						gid01 |= 1;
					if (ref_gids_perqry_arr[i][1] + 2 == ref_gids_perqry_arr[i][0])
					{
						ref_gids_perqry_arr[i][0] += 100;
						ref_gids_perqry_arr[i] = realloc(ref_gids_perqry_arr[i], ref_gids_perqry_arr[i][0] * sizeof(ref_gids_perqry_arr[i][0]));
					}
					ref_gids_perqry_arr[i][2 + ref_gids_perqry_arr[i][1]] = gid01;
					ref_gids_perqry_arr[i][1]++;
					// vector_push(&ref_gids_perqry_arr[i],&gid01);
				}
			}
			free(idx);
			qsort(ref_gids_perqry_arr[i] + 2, ref_gids_perqry_arr[i][1], sizeof(uint32_t), qsort_comparator_uint32);
			lens[i] = dedup_with_ctxobj_counts(ref_gids_perqry_arr[i] + 2, ref_gids_perqry_arr[i][1], &ctxobj_cnt_perqry_arr[i]);
		}

		for (int i = 0; i < this_block_size; i++)
		{
			int qry_gid = b * block_size + i;
			int qry_sketch_size = ref_sketch_index[qry_gid + 1] - ref_sketch_index[qry_gid];
			// sparse only code >>
			for (int l = 0; l < lens[i]; l++)
			{
				double dist = ctxobj_cnt_perqry_arr[i][l].diff_obj / ctxobj_cnt_perqry_arr[i][l].ctx_ct;
				double ani = (1 - dist);
#define VEC_GET_AS(type, vec, idx) (((type *)(vec).data)[idx])
				//			if(ctxobj_cnt_perqry_arr[i][l].ctx_ct > 20)
				//          	printf ("%d|%d|%d|%d|%d|%lf|%lf\t",qry_gid,VEC_GET_AS(uint32_t, ref_gids_perqry_arr[i], l),qry_sketch_size,ctxobj_cnt_perqry_arr[i][l].ctx_ct,ctxobj_cnt_perqry_arr[i][l].diff_obj,dist,ani);
			}
			//        printf("\n");
			//<<
		}
		//	for(int i =0 ; i< this_block_size; i++) vector_free(&ref_gids_perqry_arr[i]);//ref_gids_perqry_arr[i].size = 0;
#pragma omp parallel for num_threads(ani_opt->p) schedule(guided)
		for (int i = 0; i < this_block_size; i++)
		{
			ref_gids_perqry_arr[i][1] = 0;
			//			ref_gids_perqry_arr[i][0] = 4;
			free(ctxobj_cnt_perqry_arr[i]);
		}
		offset_gid += this_block_size;
	} // all blocks loop end

	for (int i = 0; i < block_size; i++)
		free(ref_gids_perqry_arr[i]);
	free_all(ref_gids_perqry_arr, ctxobj_cnt_perqry_arr, NULL); ////sparse only code >>
	free_all(ref_dim_sketch_stat, ref_sketch_index, qry_dim_sketch_stat, qry_sketch_index, tmp_ctxobj, NULL);
	free_read_from_file(sortedcomb_ctxgid64obj32, ctxgidobj_arr_fsize);
	fclose(fp);

	return ctxgidobj_arr_fsize;
}

// 2. inverted index(i.e. global sorted_ctxgidobj_arr) X common index(i.e. genome-wise sorted comb_sortedsketch64 ) :
//  ** the fatest method when dist matrix is sparse
void sorted_ctxgidobj_arrXcomb_sortedsketch64(unify_sketch_t *qry_result, ctxgidobj_t *ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary)
{

	uint64_t *qry_sketch_index = qry_result->sketch_index;
	uint64_t *qry_comb_sketch = qry_result->comb_sketch;
	int qry_infile_num = qry_result->infile_num;

	uint64_t gidmask = (1UL << Bitslen.gid) - 1;
	uint64_t objmask = (1UL << Bitslen.obj) - 1;
	uint32_t ref_arrlen = sort_sketch_summary->arrlen;
	int ref_infile_num = sort_sketch_summary->infile_num;
	// for self comparision only
	assert(qry_infile_num == ref_infile_num);

	uint16_t *ctx = calloc((size_t)ref_infile_num * (ref_infile_num + 1) / 2, sizeof(uint16_t));
	uint16_t *obj = calloc((size_t)ref_infile_num * (ref_infile_num + 1) / 2, sizeof(uint16_t));
#pragma omp parallel for num_threads(32) schedule(guided)
	for (uint32_t rn = 0; rn < qry_infile_num; rn++)
	{
		//		if(rn > 0) break;
		uint64_t *a = qry_comb_sketch + qry_sketch_index[rn];
		size_t a_size = qry_sketch_index[rn + 1] - qry_sketch_index[rn];
		size_t *idx = find_first_occurrences_AT_ctxgidobj_arr(a, a_size, ctxgidobj_arr, ref_arrlen);
		for (int i = 0; i < a_size; i++)
		{
			if (idx[i] == SIZE_MAX)
				continue;
			// skip conlict object;
			if ((i > 0) && ((a[i] >> Bitslen.obj) == (a[i - 1] >> Bitslen.obj)))
				continue;
			if ((i < a_size - 1) && ((a[i] >> Bitslen.obj) == (a[i + 1] >> Bitslen.obj)))
				continue;

			for (int d = idx[i];; d++)
			{
				uint32_t gid = ctxgidobj_arr[d].ctxgid & gidmask;
				if (gid > rn || (ctxgidobj_arr[d].ctxgid >> Bitslen.gid) != (a[i] >> Bitslen.obj))
					break;
				CTX(rn, gid)
				++;
				if ((a[i] & objmask) != ctxgidobj_arr[d].obj) // wrong: if(a[i] & objmask != ctxgidobj_arr[d].obj) ...
					OBJ(rn, gid)
				++;
			}
		}
		free(idx);
		/*
		#pragma omp critical
				{
					for(int i = 0 ; i <= rn;i++){
						if( CTX(rn,i) > 0 ) printf("\t%d|%d|%d|%d|%f",i,rn, CTX(rn,i),OBJ(rn,i),(float) OBJ(rn,i)/CTX(rn, i));
					}
					printf("\n");
				}
		*/
	}

	for (uint32_t rn = 0; rn < qry_infile_num; rn++)
	{
		printf("%u", rn);
		for (int i = 0; i < rn; i++)
		{
			if (CTX(rn, i) > 0)
				printf("\t%d|%d|%d|%f", i, CTX(rn, i), OBJ(rn, i), (float)OBJ(rn, i) / CTX(rn, i));
		}
		printf("\n");
	}

} // end

// 3. common index X common index (i.e. genome-wise sorted comb_sortedsketch64 ):
//     code hits: use paris wise small sortted arrays overlapping
/*  Advatages:
		** immediate output (no need precompute dist matrix )
		** most convience and memory efficience, no invert indxeing needed
	** when genome are highly similar (dense dist matrix) speed is even faster than inverted index(i.e. global sorted_ctxgidobj_arr) X common index
		** small sortted arrays overlapping is ~ 2X-3X faster than hashtable lookup based overlapping.
//...
*/


/* compute one (qn,rn); returns 1 if ani>cutoff and writes a line into ks_line */
static inline int compute_line_if_pass(    const unify_sketch_t *qry, const unify_sketch_t *ref,
    uint32_t qn, uint32_t rn, const ani_opt_t *opt,kstring_t *ks_line)
{
    uint64_t *arr_qry = qry->comb_sketch + qry->sketch_index[qn];
    size_t    len_qry = qry->sketch_index[qn + 1] - qry->sketch_index[qn];

    uint64_t *arr_ref = ref->comb_sketch + ref->sketch_index[rn];
    size_t    len_ref = ref->sketch_index[rn + 1] - ref->sketch_index[rn];

    ani_features_t f;
    get_ani_features_ctx_min_over_conflicts_a_only(arr_qry, len_qry, arr_ref, len_ref, &f);

    double af_q = (double)f.XnY_ctx / (double)len_qry;
    double af_r = (double)f.XnY_ctx / (double)len_ref;

    /* apply af cut (sets ani=0 if fail) */
    if (af_r < opt->afcut) return 0;

    double dist = get_generic_dist_from_features(&f);
    double ani  = 1.0 - dist;

    if (ani <= opt->anicut ) return 0;   /* only keep ANI > preset */

    ks_line->l = 0; /* reset single-line buffer */
    ksprintf(ks_line, "%s\t%s\t%d\t%f\t%f\t%d\t%d\t%d\t%lf\n",
             ref->gname[rn], qry->gname[qn],
             f.XnY_ctx, af_q, af_r,
             f.N_diff_obj, f.N_diff_obj_section, f.N_mut2_ctx, ani);
    return 1;
}

void comb_sortedsketch64Xcomb_sortedsketch64_filter_only(ani_opt_t *ani_opt)
{
    unify_sketch_t *qry = generic_sketch_parse(ani_opt->qrydir);
    unify_sketch_t *ref = generic_sketch_parse(ani_opt->refdir);
    if (ref->conflict)
        errx(EXIT_FAILURE, "%s(): ref '%s' contains conflicting objects!", __func__, ani_opt->refdir);

    const uint32_t Q = qry->infile_num;
    const uint32_t R = ref->infile_num;

    FILE *outfp = (ani_opt->outf[0] == '\0') ? stdout : fopen(ani_opt->outf, "w");
    if (!outfp) { perror("fopen"); return; }
    fprintf(outfp, "%s\n", print_header);

    const int P = (ani_opt->p > 0) ? ani_opt->p : 1;
    _Atomic uint32_t next_qn = 0; /* for ordered printing; comment out to allow unordered */

    #pragma omp parallel num_threads(P)
    {
        #pragma omp single nowait
        for (uint32_t qn = 0; qn < Q; ++qn) {
            #pragma omp task firstprivate(qn) shared(next_qn, qry, ref, ani_opt, outfp)
            {
                /* Per-query output buffer (one big string we’ll print once) */
                kstring_t ks_out = {0,0,0};

                if (Q >= (uint32_t)P) {
                    /* Many queries: do rn serially (best cache locality) */
                    kstring_t ks_line = {0,0,0};
                    for (uint32_t rn = 0; rn < R; ++rn) {
                        if (compute_line_if_pass(qry, ref, qn, rn, ani_opt, &ks_line))
                            kputsn(ks_line.s, ks_line.l, &ks_out);
                    }
                    free(ks_line.s);
                } else {
                    /* Few queries (e.g., Q==1): fan out rn in chunks with per-chunk local buffers */
                    #pragma omp taskgroup
                    {
                        const uint32_t CHUNK = 256; /* tune 64–1024 */
                        #pragma omp critical  /* create tasks serially */
                        for (uint32_t start = 0; start < R; start += CHUNK) {
                            const uint32_t end = (start + CHUNK < R) ? start + CHUNK : R;
                            #pragma omp task firstprivate(start,end,qn) shared(qry,ref,ani_opt,ks_out)
                            {
                                kstring_t ks_local = {0,0,0};   /* per-chunk buffer */
                                kstring_t ks_line  = {0,0,0};   /* per-line scratch */
                                for (uint32_t rn = start; rn < end; ++rn) {
                                    if (compute_line_if_pass(qry, ref, qn, rn, ani_opt, &ks_line))
                                        kputsn(ks_line.s, ks_line.l, &ks_local);
                                }
                                free(ks_line.s);

                                if (ks_local.l) {
                                    #pragma omp critical (append_qn)
                                    kputsn(ks_local.s, ks_local.l, &ks_out);
                                }
                                free(ks_local.s);
                            }
                        }
                    } /* taskgroup waits */
                }

                /* ordered print by qn (remove ticket if unordered is fine) */
                while (atomic_load_explicit(&next_qn, memory_order_acquire) != qn) {
                    #if defined(__x86_64__) || defined(__i386__)
                    __builtin_ia32_pause();
                    #endif
                }
                if (ks_out.l) fwrite(ks_out.s, 1, ks_out.l, outfp);
                atomic_fetch_add_explicit(&next_qn, 1u, memory_order_release);

                free(ks_out.s);
            } /* end task per qn */
        } /* for qn */
        /* implicit taskwait at end of parallel region */
    } /* parallel */

    if (outfp != stdout) fclose(outfp);
}

void comb_sortedsketch64Xcomb_sortedsketch64(ani_opt_t *ani_opt)
{
	unify_sketch_t *qry_result = generic_sketch_parse(ani_opt->qrydir);
	unify_sketch_t *ref_result = generic_sketch_parse(ani_opt->refdir);

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	/* choose format
		if (ani_opt->fmt)
		{ // matrix format
			for (int i = 0; i < ref_result->infile_num; i++) fprintf(outfp, "\t%s", ref_result->gname[i]);
			fprintf(outfp, "\n");
		}
		else
	*/
	fprintf(outfp, "%s\n", print_header);
	// #pragma omp parallel for num_threads(32) schedule(guided)
	for (uint32_t rn = 0; rn < ref_result->infile_num; rn++)
	{

		uint64_t *arr_ref = ref_result->comb_sketch + ref_result->sketch_index[rn];
		size_t len_ref = ref_result->sketch_index[rn + 1] - ref_result->sketch_index[rn];
		// printf(outfp,"%s", ref_result->gname[rn]);

		for (uint32_t qn = 0; qn < qry_result->infile_num; qn++)
		{
			ani_features_t ani_features;
			uint64_t *arr_qry = qry_result->comb_sketch + qry_result->sketch_index[qn];
			size_t len_qry = qry_result->sketch_index[qn + 1] - qry_result->sketch_index[qn];
			get_ani_features_from_two_sorted_ctxobj64(arr_ref, len_ref, arr_qry, len_qry, &ani_features);
			double af_qry = (double)ani_features.XnY_ctx / len_qry;
			double af_ref = (double)ani_features.XnY_ctx / len_ref;

			if (af_qry < ani_opt->afcut || af_ref < ani_opt->afcut)
				continue;
			double dist = get_generic_dist_from_features(&ani_features);
			double ani = 1 - dist;
			// double dist = get_mashD(klen, len_ref, len_qry, XnY);
			fprintf(outfp, "%s\t%s\t%d\t%f\t%f\t%d\t%d\t%d\t%lf\n", ref_result->gname[rn], qry_result->gname[qn],
					ani_features.XnY_ctx, af_qry, af_ref, ani_features.N_diff_obj, ani_features.N_diff_obj_section, ani_features.N_mut2_ctx, ani);
		}
	}
}

void check_comb_sortedsketch64(unify_sketch_t *result)
{
	for (uint32_t rn = 0; rn < result->infile_num; rn++)
	{
		uint64_t *arr = result->comb_sketch + result->sketch_index[rn];
		size_t len = result->sketch_index[rn + 1] - result->sketch_index[rn];
		for (uint32_t i = 1; i < len; i++)
		{
			if (arr[i] < arr[i - 1])
				err(EXIT_FAILURE, "%s(): %dth genome %dth kmer < %dth kmer (%lx<%lx)", __func__, rn, i, i - 1, arr[i], arr[i - 1]);
		}
	}
}

//

/**
 * Finds the first occurrence index in large sorted array `b` for each element in a small sorted array `a`.
 *
 * @param a       Sorted array of uint64_t elements (ascending order)
 * @param a_size  Number of elements in array `a`
 * @param b       Sorted array of uint64_t elements (ascending order)
 * @param b_size  Number of elements in array `b`
 * @return        Array of indices (size_t*) where indices[i] = first occurrence of a[i] in b,
 *                SIZE_MAX if not found. Caller must free() the returned array.
 */
size_t *find_first_occurrences_AT_ctxgidobj_arr(const uint64_t *a, size_t a_size,
												const ctxgidobj_t *b, size_t b_size)
{
	size_t *indices = malloc(a_size * sizeof(size_t));
	if (!indices)
		return NULL;
	size_t low = 0; // Track lower bound for binary search
	int nobjbits = Bitslen.obj;

	for (size_t i = 0; i < a_size; ++i)
	{
		const uint64_t target = a[i] >> nobjbits; // (2*(2*holen + iolen ));
		// when conflict objects are kept, skip searching if target[i+1] == target[i].
		if (i > 0 && target == a[i - 1] >> nobjbits)
		{
			indices[i] = indices[i - 1]; // Use the previous index if the current target is the same
			continue;
		}

		size_t high = b_size;
		// Leftmost binary search within [low, high)
		while (low < high)
		{
			size_t mid = low + (high - low) / 2;
			if ((b[mid].ctxgid >> GID_NBITS) < target)
			{
				low = mid + 1;
			}
			else
			{
				high = mid;
			}
		}

		// Check if target was found: modified from: if (low < b_size && b[low] == target) {
		if (low < b_size && (b[low].ctxgid >> GID_NBITS == target))
			indices[i] = low;
		else
			indices[i] = SIZE_MAX; // Not found
	}

	return indices;
}

void ani_block_print(
	int ref_infile_num, int qry_gid_offset, int this_block_size,
	uint64_t *ref_sketch_index, uint64_t *qry_sketch_index,
	ctx_mut2_t *ctx, obj_section_t *obj,
	char (*refname)[PATHLEN], char (*qryfname)[PATHLEN],
	uint32_t *num_passid_block, idani_t **sort_idani_block,
	FILE *outfp, ani_opt_t *ani_opt, int matrix_mode)
{
	ani_features_t ani_features;

	for (int i = 0; i < this_block_size; i++)
	{
		int qry_gid = qry_gid_offset + i;
		int qry_sketch_size = qry_sketch_index[qry_gid + 1] - qry_sketch_index[qry_gid];

		if (matrix_mode)
		{
			fprintf(outfp, "%s", qryfname[qry_gid]);
		}

		int loop_j = matrix_mode ? ref_infile_num : num_passid_block[i];
		for (int n = 0; n < loop_j; n++)
		{
			int j = matrix_mode ? n : sort_idani_block[i][n].id;
			ani_features.XnY_ctx = MCTX(ref_infile_num, i, j).num_ctx;
			ani_features.N_diff_obj_section = MOBJ(ref_infile_num, i, j).diff_obj_section;
			ani_features.N_mut2_ctx = MCTX(ref_infile_num, i, j).num_mut2_ctx;
			ani_features.N_diff_obj = MOBJ(ref_infile_num, i, j).diff_obj;

			int ref_sketch_size = ref_sketch_index[j + 1] - ref_sketch_index[j];
			float af_qry = (float)ani_features.XnY_ctx / qry_sketch_size;
			float af_ref = (float)ani_features.XnY_ctx / ref_sketch_size;
			float min_af = af_qry < af_ref ? af_qry : af_ref;

			select_metrics_dist[0] = min_af < ani_opt->afcut ? ani_opt->e : get_generic_dist_from_features(&ani_features);
			select_metrics_dist[1] = get_mashD(Bitslen.ctx / 2, ref_sketch_size, qry_sketch_size, ani_features.XnY_ctx);
			select_metrics_dist[2] = get_aafD(Bitslen.ctx / 2, ref_sketch_size, qry_sketch_size, ani_features.XnY_ctx);
			select_metrics_dist[3] = min_af < ani_opt->afcut ? select_metrics_dist[1] : select_metrics_dist[0];
			select_metrics_dist[4] = min_af < ani_opt->afcut ? select_metrics_dist[2] : select_metrics_dist[0];

			double dist = select_metrics_dist[abs(ani_opt->s) - 1];
			double ani = 1 - dist;
			double metric = ani_opt->s > 0 ? ani : dist;

			if (matrix_mode)
			{
				fprintf(outfp, "\t%lf", metric);
			}
			else
			{
				fprintf(outfp, "%s\t%s\t%d\t%f\t%f\t%d\t%d\t%d\t%lf\t%lf\n",
						qryfname[qry_gid], refname[j], ani_features.XnY_ctx, af_qry, af_ref,
						MOBJ(ref_infile_num, i, j).diff_obj,
						MOBJ(ref_infile_num, i, j).diff_obj_section,
						MCTX(ref_infile_num, i, j).num_mut2_ctx, ani, dist);
			}
		}
		if (matrix_mode)
		{
			fprintf(outfp, "\n");
		}
	}
}

void simple_sortedsketch64Xcomb_sortedsketch64(simple_sketch_t *simple_sketch, infile_tab_t *genomes_infiletab, ani_opt_t *ani_opt)
{
	uint64_t *arr_ref = simple_sketch->comb_sketch;
	if (arr_ref == NULL)
		err(EXIT_FAILURE, "%s(): simple_sketch->comb_sketch is NULL", __func__);

	if (simple_sketch->infile_num != genomes_infiletab->infile_num)
		err(EXIT_FAILURE, "%s(): infile_num mismatch: %d vs %d", __func__, simple_sketch->infile_num, genomes_infiletab->infile_num);
	size_t ref_sketch_size = simple_sketch->sketch_index[1] - simple_sketch->sketch_index[0];
	assert(simple_sketch->infile_num == genomes_infiletab->infile_num);

	for (uint32_t qn = 1; qn < simple_sketch->infile_num; qn++)
	{
		ani_features_t ani_features;
		uint64_t *arr_qry = simple_sketch->comb_sketch + simple_sketch->sketch_index[qn];
		size_t qry_sketch_size = simple_sketch->sketch_index[qn + 1] - simple_sketch->sketch_index[qn];

		get_ani_features_from_two_sorted_ctxobj64(arr_ref, ref_sketch_size, arr_qry, qry_sketch_size, &ani_features);

		double af_qry = (double)ani_features.XnY_ctx / qry_sketch_size;

		double af_ref = (double)ani_features.XnY_ctx / ref_sketch_size;

		float min_af = af_qry < af_ref ? af_qry : af_ref;

		select_metrics_dist[0] = min_af < ani_opt->afcut ? ani_opt->e : get_generic_dist_from_features(&ani_features);
		select_metrics_dist[1] = get_mashD(Bitslen.ctx / 2, ref_sketch_size, qry_sketch_size, ani_features.XnY_ctx);
		select_metrics_dist[2] = get_aafD(Bitslen.ctx / 2, ref_sketch_size, qry_sketch_size, ani_features.XnY_ctx);
		select_metrics_dist[3] = min_af < ani_opt->afcut ? select_metrics_dist[1] : select_metrics_dist[0];
		select_metrics_dist[4] = min_af < ani_opt->afcut ? select_metrics_dist[2] : select_metrics_dist[0];

		double dist = select_metrics_dist[abs(ani_opt->s) - 1];
		double ani = 1 - dist;
		printf("%s\t%s\t%d\t%f\t%f\t%d\t%d\t%d\t%lf\n", genomes_infiletab->organized_infile_tab[0].fpath, genomes_infiletab->organized_infile_tab[qn].fpath,
			   ani_features.XnY_ctx, af_qry, af_ref, ani_features.N_diff_obj, ani_features.N_diff_obj_section, ani_features.N_mut2_ctx, ani);
	}
}

// fencepost search method
#include <limits.h>

static inline uint32_t top_k_bits_u64(uint64_t x, int k)
{
	if (k <= 0)
		return 0u;
	if (k >= 64)
		return (uint32_t)x; /* defensive; we clamp k well below 64 */
	return (uint32_t)(x >> (64 - k));
}

int kssd_choose_k_fenceposts(size_t b_size, size_t a_size)
{
	if (a_size == 0 || b_size == 0)
		return 14;
	double ratio = (double)b_size / (double)a_size;			   /* ≈ m = b/a */
	int k = (ratio > 1.0) ? (int)floor(log2(ratio) + 0.5) : 0; /* ~log2(m) */
	if (k < 12)
		k = 12; /* 4K buckets  (~32 KB of fenceposts) */
	if (k > 20)
		k = 20; /* 1M buckets  (~8  MB of fenceposts, 64-bit size_t) */
	return k;
}

int kssd_build_fenceposts_ctxgid(const ctxgidobj_t *b, size_t b_size, int k, size_t *F)
{
	if (!b || !F || k < 0 || k > 32)
		return -1;
	const size_t buckets = (size_t)1u << k;

	for (size_t t = 0; t <= buckets; ++t)
		F[t] = b_size;

	size_t next = 0;
	for (size_t i = 0; i < b_size; ++i)
	{
		uint64_t key = (uint64_t)b[i].ctxgid >> GID_NBITS; /* effective key */
		uint32_t topk = top_k_bits_u64(key, k);
		while (next <= topk)
			F[next++] = i;
		if (next > buckets)
			break;
	}
	while (next <= buckets)
		F[next++] = b_size;
	return 0;
}

static inline size_t lb_in_bucket_ctxgid(const ctxgidobj_t *b, const size_t *F, int k,
										 uint64_t target_key)
{
	const size_t buckets = (size_t)1u << k;
	uint32_t t = top_k_bits_u64(target_key, k);
	if (t > buckets)
		t = (uint32_t)buckets; /* defensive */

	size_t lo = F[t];
	size_t hi = F[t + 1];

	while (lo < hi)
	{
		size_t mid = lo + ((hi - lo) >> 1);
		uint64_t mid_key = (uint64_t)b[mid].ctxgid >> GID_NBITS;
		if (mid_key < target_key)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo; /* first position with key >= target_key (may be hi) */
}

size_t *kssd_find_first_occurrences_fenceposts(const uint64_t *a, size_t a_size,
											   const ctxgidobj_t *b, size_t b_size,
											   const size_t *F, int k, unsigned nobjbits)
{
	if (!a || !b || !F)
		return NULL;

	size_t *idx = (size_t *)malloc(a_size * sizeof *idx);
	if (!idx)
		return NULL;

	uint64_t prev_key = UINT64_MAX;
	size_t prev_idx = SIZE_MAX;

	for (size_t i = 0; i < a_size; ++i)
	{
		uint64_t target_key = a[i] >> nobjbits;

		if (i && target_key == prev_key)
		{ /* reuse for duplicates in a */
			idx[i] = prev_idx;
			continue;
		}

		size_t pos = lb_in_bucket_ctxgid(b, F, k, target_key);

		if (pos < b_size && (((uint64_t)b[pos].ctxgid >> GID_NBITS) == target_key))
			idx[i] = pos;
		else
			idx[i] = SIZE_MAX;

		prev_key = target_key;
		prev_idx = idx[i];
	}
	return idx;
}
