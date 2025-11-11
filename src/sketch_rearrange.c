#include "sketch_rearrange.h"

uint64_t ctxmask, tupmask, ho_mask_len, hc_mask_len, io_mask_len, ho_mask_left, hc_mask_left, io_mask, hc_mask_right, ho_mask_right;
uint8_t iolen, klen, hclen, holen;
uint32_t gid_mask;
bitslen_t Bitslen;
void const_comask_init(dim_sketch_stat_t *lco_stat_val)
{
    //  init all public vars ;
    klen = lco_stat_val->klen;
    if (lco_stat_val->coden_len > 0)
    {
        Bitslen.ctx = 4 * lco_stat_val->coden_len;
        ctxmask = generate_coden_pattern64();
    }
    else
    {
        holen = lco_stat_val->holen;
        hclen = lco_stat_val->hclen;
        iolen = klen - 2 * (hclen + holen);

        Bitslen.ctx = 4 * hclen;
        ho_mask_len = holen == 0 ? 0 : UINT64_MAX >> (64 - 2 * holen);
        hc_mask_len = hclen == 0 ? 0 : UINT64_MAX >> (64 - 2 * hclen);
        io_mask_len = iolen == 0 ? 0 : UINT64_MAX >> (64 - 2 * iolen); // UINT64_MAX >> 64 is undefined  not 0

        ho_mask_left = ho_mask_len << (2 * (klen - holen));          //(2*(holen + hclen + iolen + hclen));
        hc_mask_left = hc_mask_len << (2 * (holen + hclen + iolen)); // 2*(holen + hclen + iolen));
        io_mask = io_mask_len << (2 * (holen + hclen));
        hc_mask_right = hc_mask_len << (2 * holen);
        ho_mask_right = ho_mask_len;

        ctxmask = hc_mask_left | hc_mask_right;
    }
    Bitslen.obj = 2 * klen - Bitslen.ctx;
    Bitslen.gid = GID_NBITS;
    tupmask = klen == 0 ? 0 : UINT64_MAX >> (64 - 2 * klen);
    gid_mask = (1U << GID_NBITS) - 1;
}

void sketch64_2ctxobj64(uint64_t *sketch64, uint32_t arrlen)
{
#pragma omp parallel for num_threads(32) schedule(guided)
    for (uint32_t ri = 0; ri < arrlen; ri++)
    {
        sketch64[ri] = uint64kmer2generic_ctxobj(sketch64[ri]);
    }
}

uint96_t *sketch64_2uint96co(uint64_t *sketch_index, uint64_t *sketch64, int infile_num, uint32_t arrlen)
{
    //  int iolen = klen - 2*(hclen + holen);
    uint96_t *uint96co = malloc(arrlen * sizeof(uint96_t));
#pragma omp parallel for num_threads(32) schedule(guided)
    for (uint32_t rn = 0; rn < infile_num; rn++)
    {
        for (uint32_t ri = sketch_index[rn]; ri < sketch_index[rn + 1]; ri++)
        {
            uint96co[ri] = uint64_kmer2uint96(sketch64[ri], rn);
        }
    }
    return uint96co;
}

ctxgidobj_t *ctxobj64_2ctxgidobj(uint64_t *sketch_index, uint64_t *ctxobj64, int infile_num, uint32_t arrlen)
{
    ctxgidobj_t *ctxgidobj = malloc(arrlen * sizeof(ctxgidobj_t));
    uint64_t obj_len_mask = (1LU << Bitslen.obj) - 1;
#pragma omp parallel for num_threads(32) schedule(guided)
    for (uint32_t rn = 0; rn < infile_num; rn++)
    {
        for (uint32_t ri = sketch_index[rn]; ri < sketch_index[rn + 1]; ri++)
        {
            ctxgidobj[ri] = uint64_ctxobj2ctxgidobj96(ctxobj64[ri], rn, obj_len_mask);
        }
    }
    return ctxgidobj;
}

uint64kmer2generic_ctxobj_fn uint64kmer2generic_ctxobj = NULL;
void set_uint64kmer2generic_ctxobj(bool is_coden_ctxobj_pattern)
{
    
    uint64kmer2generic_ctxobj = is_coden_ctxobj_pattern ? reorder_unituple_by_coden_pattern64 : uint64_kmer2ctxobj;
}
