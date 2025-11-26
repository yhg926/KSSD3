#define _GNU_SOURCE   // must come before any system header
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "global_basic.h"
#include "command_sketch.h"
#include "sketch_rearrange.h"
#include "kssdlib_sort.h"
#include <time.h>

#ifndef likely
#define likely(x) __builtin_expect(!!(x), 1)
#endif
#ifndef unlikely
#define unlikely(x) __builtin_expect(!!(x), 0)
#endif

// shared global vars
const char sketch_stat[] = "lcofiles.stat";
const char sketch_suffix[] = "lco"; // long co, 64bits
const char combined_sketch_suffix[] = "comblco";
const char idx_sketch_suffix[] = "comblco.index";
const char combined_ab_suffix[] = "comblco.a";
const char sorted_comb_ctxgid64obj32[] = "sortedcomb_ctxgid64obj32";
// public vars shared across files
uint32_t FILTER, hash_id;
dim_sketch_stat_t comblco_stat_one, comblco_stat_it;
// tmp container in this scope
static size_t file_size;
static char tmp_fname[PATHLEN + 20];
static struct stat tmpstat;

// functions
void test_read_genomes2mem2sortedctxobj64(sketch_opt_t *sketch_opt_val,infile_tab_t *infile_stat,int batch_size);
uint32_t get_sketching_id(uint32_t hclen, uint32_t holen, uint32_t iolen, uint32_t drfold, uint32_t FILTER)
{
    return GET_SKETCHING_ID(hclen, holen, iolen, drfold, FILTER);
}

void compute_sketch(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat)
{
    if (sketch_opt_val->split_mfa)
    { // mfa files parse
        mfa2sortedctxobj64_v2(sketch_opt_val, infile_stat);
        return;
    }
    // test_read_genomes2mem2sortedctxobj64(sketch_opt_val, infile_stat, 1000);
    // read_genomes2mem2sortedctxobj64(sketch_opt_val, infile_stat, 1000);

    // Decide mode based on #files vs. threads
    if (infile_stat->infile_num >= sketch_opt_val->p)         // Mode-A: per-file parallel
        sketch_many_files_in_parallel(sketch_opt_val, infile_stat, 1024);
    else         // Mode-B: in-file parallel
        sketch_few_files_with_intrafile_parallel(sketch_opt_val, infile_stat, 4096);

    return;
    // 20150724: i am not sure if below code is still needed, since we have read_genomes2mem2sortedctxobj64() to handle all cases
    // normal sketching mode
    uint64_t *tmp_ct_list = calloc(infile_stat->infile_num + 1, sizeof(uint64_t));
#pragma omp parallel for num_threads(sketch_opt_val->p) schedule(guided)
    for (int i = 0; i < infile_stat->infile_num; i++)
    {
        // sketching all genomes individually
        tmp_ct_list[i + 1] = seq2ht_sortedctxobj64((infile_stat->organized_infile_tab)[i].fpath, format_string("%s/%d.%s", sketch_opt_val->outdir, i, sketch_suffix), sketch_opt_val->abundance, sketch_opt_val->kmerocrs);
        /* may consider other sketching methods, i.e. sorting arr and dedup/count instead of hashtable :
         *  	reads2sketch64(...)
         *		 seq2sortedsketch64()
         *   	opt_seq2sortedsketch64();
         */
        printf("\r%dth/%d genome sketching %s completed!\tsketch size=%lu", i + 1, infile_stat->infile_num,
               (infile_stat->organized_infile_tab)[i].fpath, tmp_ct_list[i + 1]);
    }
    printf("\n");
    // combine *.lco to comblco
    combine_lco(sketch_opt_val, infile_stat);
    // write comblco.index
    for (int i = 0; i < infile_stat->infile_num; i++)
        tmp_ct_list[i + 1] += tmp_ct_list[i];
    write_to_file(test_create_fullpath(sketch_opt_val->outdir, idx_sketch_suffix), tmp_ct_list, sizeof(tmp_ct_list[0]) * (infile_stat->infile_num + 1));
    free(tmp_ct_list);
    // write stat file
    write_sketch_stat(sketch_opt_val->outdir, infile_stat);
};

void combine_lco(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat)
{
    // combine *.lco to comblco
    int i = 0;
    char *comb_ab_fn, *comb_sketch_fn;
    FILE *comb_sketch_fp, *comb_ab_fp;
    comb_sketch_fn = format_string("%s/%s", sketch_opt_val->outdir, combined_sketch_suffix);
    sprintf(tmp_fname, "%s/%d.%s", sketch_opt_val->outdir, i, sketch_suffix);
    if (rename(tmp_fname, comb_sketch_fn))
        err(errno, "%s():%s rename error", __func__, tmp_fname);
    if ((comb_sketch_fp = fopen(comb_sketch_fn, "ab")) == NULL)
        err(errno, "%s() open file error: %s", __func__, comb_sketch_fn);
    if (sketch_opt_val->abundance)
    {
        comb_ab_fn = format_string("%s.a", comb_sketch_fn);
        sprintf(tmp_fname, "%s/%d.%s.a", sketch_opt_val->outdir, i, sketch_suffix);
        if (rename(tmp_fname, comb_ab_fn))
            err(errno, "%s():%s rename error", __func__, tmp_fname);
        if ((comb_ab_fp = fopen(comb_ab_fn, "ab")) == NULL)
            err(errno, "%s() open file error: %s", __func__, comb_ab_fn);
    }

    for (int i = 1; i < infile_stat->infile_num; i++)
    {
        sprintf(tmp_fname, "%s/%d.%s", sketch_opt_val->outdir, i, sketch_suffix);
        uint64_t *mem_lco = read_from_file(tmp_fname, &file_size);
        fwrite(mem_lco, file_size, 1, comb_sketch_fp);
        remove(tmp_fname);
        free(mem_lco);
        // abundance
        if (!sketch_opt_val->abundance)
            continue;
        sprintf(tmp_fname, "%s/%d.%s.a", sketch_opt_val->outdir, i, sketch_suffix);
        uint32_t *mem_ab = read_from_file(tmp_fname, &file_size);
        fwrite(mem_ab, file_size, 1, comb_ab_fp);
        remove(tmp_fname);
        free(mem_ab);
    }
    fclose(comb_sketch_fp);
    if (sketch_opt_val->abundance)
        fclose(comb_ab_fp);
}

void gen_inverted_index4comblco(const char *refdir)
{

    unify_sketch_t *ref_result = generic_sketch_parse(refdir);
    const_comask_init(&ref_result->stats.lco_stat_val);

    uint64_t sketch_size = ref_result->sketch_index[ref_result->infile_num];
    if (sketch_size >= UINT32_MAX)
        err(EXIT_FAILURE, "%s():sketch_index maximun %lu exceed UINT32_MAX %u", __func__, sketch_size, UINT32_MAX);
    if (ref_result->infile_num >= (1 << GID_NBITS))
        err(EXIT_FAILURE, "%s(): genome numer %d exceed maximum:%u", __func__, ref_result->infile_num, 1 << GID_NBITS);
    if (GID_NBITS + 4 * hclen > 64)
        err(EXIT_FAILURE, "%s(): context_bits_len(%d)+gid_bits_len(%d) exceed 64", __func__, 4 * hclen, GID_NBITS);
    ctxgidobj_t *ctxgidobj = ctxobj64_2ctxgidobj(ref_result->sketch_index, ref_result->comb_sketch, ref_result->infile_num, sketch_size);
    free_unify_sketch(ref_result);
    ctxgidobj_sort_array(ctxgidobj, sketch_size);
    // printf("sketch_size=%lu\t%d\t%d\n",sizeof(ctxgidobj[0]),sizeof(ctxgidobj_t),sketch_size);
    write_to_file(format_string("%s/%s", refdir, sorted_comb_ctxgid64obj32), ctxgidobj, sizeof(ctxgidobj[0]) * sketch_size);
    free(ctxgidobj);
}

void write_sketch_stat(const char *outdir, infile_tab_t *infile_stat)
{
    char (*tmpname)[PATHLEN] = malloc(infile_stat->infile_num * PATHLEN);
    for (int i = 0; i < infile_stat->infile_num; i++)
        memcpy(tmpname[i], (infile_stat->organized_infile_tab)[i].fpath, PATHLEN);
    concat_and_write_to_file(test_create_fullpath(outdir, sketch_stat), &comblco_stat_one, sizeof(comblco_stat_one), tmpname, infile_stat->infile_num * PATHLEN);
    free(tmpname);
}

KSEQ_INIT(gzFile, gzread)
// keep k-mer using hashtable is faster(~16s) than using vector and sort(~24) per 1k genomes, probabaly due to cache competition
int seq2ht_sortedctxobj64(char *seqfname, char *outfname, bool abundance, int n)
{
    int TL = klen;
    gzFile infile = gzopen(seqfname, "r");
    if (!infile)
        err(errno, "reads2sketch64(): Cannot open file %s", seqfname);
    kseq_t *seq = kseq_init(infile);
    khash_t(sort64) *h = kh_init(sort64);
    uint64_t tuple, crvstuple, unituple, basenum, unictx;
    uint32_t len_mv = 2 * TL - 2;
    uint32_t sketch_size = 0;

    while (kseq_read(seq) >= 0)
    {
        const char *s = seq->seq.s;
        if (seq->seq.l < TL)
            continue;
        int base = 0;

        for (int pos = 0; pos < seq->seq.l; pos++)
        { // for(int pos = TL; pos < seq->seq.l ; pos++)
            if (Basemap[(unsigned short)s[pos]] == DEFAULT)
            {
                base = 0;
                continue;
            }
            basenum = Basemap[(unsigned short)s[pos]];
            tuple = ((tuple << 2) | basenum);
            crvstuple = ((crvstuple >> 2) | ((basenum ^ 3LLU) << len_mv));
            if (++base < TL)
                continue;
            // if base >=TL, namely, contiue ACGT TL-mer
            unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
            unictx = unituple & ctxmask;
            if (SKETCH_HASH(unictx) > FILTER)
                continue;
            int ret;
            khint_t key = kh_put(sort64, h, uint64kmer2generic_ctxobj(unituple), &ret);

            if (ret)
            {
                kh_value(h, key) = 1;
                sketch_size++;
            }
            else
            {
                kh_value(h, key)++;
            }
        } // for line
    }; // while
    kseq_destroy(seq);
    gzclose(infile);
    // get sketch sorted
    SortedKV_Arrays_t lco_ab = sort_khash_u64(h);
    //	for(int i = 0; i<lco_ab.len;i++){		printf("%d\t%d\t%lx\n",i,lco_ab.len,lco_ab.keys[i]);}
    kh_destroy(sort64, h);
    if (n > 1)
        filter_n_SortedKV_Arrays(&lco_ab, n);
    write_to_file(outfname, lco_ab.keys, lco_ab.len * sizeof(lco_ab.keys[0]));
    if (abundance)
        write_to_file(format_string("%s.a", outfname), lco_ab.values, lco_ab.len * sizeof(lco_ab.values[0]));
    free_all(lco_ab.keys, lco_ab.values, NULL);
    return lco_ab.len; // sketch_size;
}
//
int seq2sortedsketch64(char *seqfname, char *outfname, bool abundance, int n)
{

    gzFile infile = gzopen(seqfname, "r");
    if (!infile)
        err(errno, "reads2sketch64(): Cannot open file %s", seqfname);
    kseq_t *seq = kseq_init(infile);
    Vector raw_sketch;
    vector_init(&raw_sketch, sizeof(uint64_t));

    uint64_t tuple, crvstuple, unituple, basenum, unictx;
    uint32_t len_mv = 2 * klen - 2;

    while (kseq_read(seq) >= 0)
    {
        const char *s = seq->seq.s;

        if (seq->seq.l < klen)
            continue;
        int base = 0;

        for (int pos = 0; pos < seq->seq.l; pos++)
        { // for(int pos = klen; pos < seq->seq.l ; pos++)
            if (Basemap[(unsigned short)s[pos]] == DEFAULT)
            {
                base = 0;
                continue;
            }
            basenum = Basemap[(unsigned short)s[pos]];
            tuple = ((tuple << 2) | basenum);
            crvstuple = ((crvstuple >> 2) | ((basenum ^ 3LLU) << len_mv));
            if (++base < klen)
                continue;
            // if base >=klen, namely, contiue ACGT klen-mer
            unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
            unictx = unituple & ctxmask;
            if (SKETCH_HASH(unictx) > FILTER)
                continue;

            unituple = uint64kmer2generic_ctxobj(unituple);
            vector_push(&raw_sketch, &unituple);

        } // for line
    }; // while
    gzclose(infile);

    // write sketch and abundance
    uint64_t *mem_lco = raw_sketch.data;
    qsort(mem_lco, raw_sketch.size, sizeof(uint64_t), qsort_comparator_uint64);
    uint32_t *mem_ab; // if (abundance) mem_ab = malloc(sketch_size * sizeof(uint32_t));
    uint32_t sketch_size, kmer_ct;
    sketch_size = kmer_ct = dedup_with_counts(mem_lco, raw_sketch.size, &mem_ab);

    if (n > 1)
    {
        kmer_ct = 0;
        for (size_t i = 0; i < sketch_size; i++)
        {
            if (mem_ab[i] >= n)
            {
                mem_lco[kmer_ct] = mem_lco[i]; // Shift unique element forward
                mem_ab[kmer_ct] = mem_ab[i];
                kmer_ct++;
            }
        }
    }

    write_to_file(outfname, mem_lco, kmer_ct * sizeof(uint64_t));
    vector_free(&raw_sketch);

    if (abundance)
    {
        write_to_file(format_string("%s.a", outfname), mem_ab, kmer_ct * sizeof(uint32_t));
        free(mem_ab);
    }
    return kmer_ct; // sketch_size;
}


int merge_comblco(sketch_opt_t *sketch_opt_val)
{

    void *mem_stat = read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[0], sketch_stat), &file_size);
    memcpy(&comblco_stat_one, mem_stat, sizeof(comblco_stat_one));
    comblco_stat_one.infile_num = 0;
    char (*tmpname)[PATHLEN] = malloc(PATHLEN);
    FILE *merge_out_fp = fopen(test_create_fullpath(sketch_opt_val->outdir, combined_sketch_suffix), "wb");
    if (merge_out_fp == NULL)
        err(errno, "%s():%s/%s", __func__, sketch_opt_val->outdir, combined_sketch_suffix);

    uint64_t *index_arry = malloc(sizeof(uint64_t) * 100);
    index_arry[0] = 0;
    for (int i = 0; i < sketch_opt_val->num_remaining_args; i++)
    {
        // read stat
        void *mem_stat_it = read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[i], sketch_stat), &file_size);
        memcpy(&comblco_stat_it, mem_stat_it, sizeof(comblco_stat_it));
        if (comblco_stat_it.hash_id != comblco_stat_one.hash_id)
            err(EINVAL, "%uth %s hashid: %u != %u ", i, sketch_opt_val->remaining_args[i], comblco_stat_it.hash_id, comblco_stat_one.hash_id);
        if (comblco_stat_it.koc == 0)
            comblco_stat_one.koc = 0;
        // collect filenames
        tmpname = realloc(tmpname, PATHLEN * (comblco_stat_one.infile_num + comblco_stat_it.infile_num));
        memcpy(tmpname + comblco_stat_one.infile_num, mem_stat_it + sizeof(comblco_stat_it), PATHLEN * comblco_stat_it.infile_num);
        // set index
        uint64_t *mem_index_it = read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[i], idx_sketch_suffix), &file_size);
        index_arry = (uint64_t *)realloc(index_arry, sizeof(uint64_t) * (comblco_stat_one.infile_num + comblco_stat_it.infile_num + 1));
        for (int j = 1; j < comblco_stat_it.infile_num + 1; j++)
            index_arry[comblco_stat_one.infile_num + j] = index_arry[comblco_stat_one.infile_num] + mem_index_it[j];
        // add file num
        comblco_stat_one.infile_num += comblco_stat_it.infile_num;
        // write combined_sketch_suffix
        uint64_t *mem_comblco = read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[i], combined_sketch_suffix), &file_size);
        fwrite(mem_comblco, file_size, 1, merge_out_fp);
        free_all(mem_stat_it, mem_index_it, mem_comblco, NULL);
    }
    fclose(merge_out_fp); // write combined_sketch_suffix complete

    if (comblco_stat_it.koc)
    {
        merge_out_fp = fopen(test_create_fullpath(sketch_opt_val->outdir, combined_ab_suffix), "wb");
        if (merge_out_fp == NULL)
            err(errno, "%s():%s/%s", __func__, sketch_opt_val->outdir, combined_ab_suffix);
        for (int i = 0; i < sketch_opt_val->num_remaining_args; i++)
        {
            uint64_t *mem_ab = read_from_file(test_get_fullpath(sketch_opt_val->remaining_args[i], combined_ab_suffix), &file_size);
            fwrite(mem_ab, file_size, 1, merge_out_fp);
            free(mem_ab);
        }
        fclose(merge_out_fp);
    }
    // write index and stat file
    write_to_file(test_create_fullpath(sketch_opt_val->outdir, idx_sketch_suffix), index_arry, sizeof(index_arry[0]) * (comblco_stat_one.infile_num + 1));
    concat_and_write_to_file(test_create_fullpath(sketch_opt_val->outdir, sketch_stat), &comblco_stat_one, sizeof(comblco_stat_one), tmpname, PATHLEN * (comblco_stat_one.infile_num));

    return comblco_stat_one.infile_num;
}



KHASH_MAP_INIT_INT64(kmer_hash, int)
void mfa2sortedctxobj64(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat)
{

    uint64_t tuple, crvstuple, unituple, basenum, unictx, totle_sketch_size = 0;
    uint32_t len_mv = 2 * klen - 2; // uint32_t saved_genome_num = infile_stat->infile_num, genome_num = 0;
    Vector sketch_index;
    vector_init(&sketch_index, sizeof(uint64_t));
    vector_push(&sketch_index, &totle_sketch_size);
    int tmpname_size_alloc = sketch_index.capacity;
    char (*tmpname)[PATHLEN] = malloc(tmpname_size_alloc * PATHLEN);
    if (!tmpname)
        err(errno, "%s():Memory allocation failed for tmpname", __func__);
    FILE *comb_sketch_fp;
    if ((comb_sketch_fp = fopen(format_string("%s/%s", sketch_opt_val->outdir, combined_sketch_suffix), "wb")) == NULL)
        err(errno, "%s() open file error: %s/%s", __func__, sketch_opt_val->outdir, combined_sketch_suffix);

    for (int i = 0; i < infile_stat->infile_num; i++)
    {
        char *seqfname = (infile_stat->organized_infile_tab)[i].fpath;
        gzFile infile = gzopen(seqfname, "r");
        if (!infile)
            err(errno, "reads2sketch64(): Cannot open file %s", seqfname);
        kseq_t *seq = kseq_init(infile);
        while (kseq_read(seq) >= 0 && seq->seq.l > klen)
        {
            // seq name handle
            if (sketch_index.size >= tmpname_size_alloc)
            {
                tmpname_size_alloc += 1000;
                tmpname = realloc(tmpname, tmpname_size_alloc * PATHLEN);
            }
            replace_special_chars_with_underscore(seq->name.s);
            strncpy(tmpname[sketch_index.size - 1], seq->name.s, PATHLEN);
            // sketching each fasta sequence
            khash_t(sort64) *h = kh_init(sort64);
            const char *s = seq->seq.s;
            int base = 0;
            uint32_t sketch_size = 0;

            for (int pos = 0; pos < seq->seq.l; pos++)
            {
                if (Basemap[(unsigned short)s[pos]] == DEFAULT)
                {
                    base = 0;
                    continue;
                }
                basenum = Basemap[(unsigned short)s[pos]];
                tuple = ((tuple << 2) | basenum);
                crvstuple = ((crvstuple >> 2) | ((basenum ^ 3LLU) << len_mv));
                if (++base < klen)
                    continue;

                unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
                unictx = unituple & ctxmask;
                if (SKETCH_HASH(unictx) > FILTER)
                    continue;
                int ret;
                khint_t key = kh_put(sort64, h, uint64kmer2generic_ctxobj(unituple), &ret);

                if (ret)
                    kh_value(h, key) = 1;

            } // for line
            SortedKV_Arrays_t lco_ab = sort_khash_u64(h);
            totle_sketch_size += lco_ab.len;
            vector_push(&sketch_index, &totle_sketch_size);
            kh_destroy(sort64, h);
            fwrite(lco_ab.keys, sizeof(lco_ab.keys[0]), lco_ab.len, comb_sketch_fp);
            free_all(lco_ab.keys, lco_ab.values, NULL);
        } // while
        kseq_destroy(seq);
        gzclose(infile);
        printf("\r%dth/%d multifasta sketching %s completed!\t #genomes=%lu", i + 1,
               infile_stat->infile_num, (infile_stat->organized_infile_tab)[i].fpath, sketch_index.size - 1);
        if (i == infile_stat->infile_num - 1)
            printf("\n");
    } // for infile
    fclose(comb_sketch_fp);
    // write index	and stat file
    write_to_file(test_create_fullpath(sketch_opt_val->outdir, idx_sketch_suffix), sketch_index.data, sizeof(uint64_t) * sketch_index.size);
    comblco_stat_one.infile_num = sketch_index.size - 1;
    concat_and_write_to_file(test_create_fullpath(sketch_opt_val->outdir, sketch_stat),
                             &comblco_stat_one, sizeof(comblco_stat_one), tmpname, comblco_stat_one.infile_num * PATHLEN);

    vector_free(&sketch_index);
    free(tmpname);
} // mfa2sortedctxobj64()

void read_genomes2mem2sortedctxobj64(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat, int batch_size)
{
    //   printf("%lx\n", ctxmask);
    uint32_t len_mv = 2 * klen - 2;
    uint64_t *sketch_index = calloc(infile_stat->infile_num + 1, sizeof(uint64_t));
    int *gseq_nums = calloc(batch_size + 1, sizeof(int));
    uint64_t **batch_sketches = malloc(batch_size * sizeof(uint64_t *));
    Vector all_reads;
    vector_init(&all_reads, sizeof(char *));

    FILE *comb_sketch_fp;
    if ((comb_sketch_fp = fopen(format_string("%s/%s", sketch_opt_val->outdir, combined_sketch_suffix), "wb")) == NULL)
        err(errno, "%s() open file error: %s/%s", __func__, sketch_opt_val->outdir, combined_sketch_suffix);
    for (int infile_num_p = 0; infile_num_p < infile_stat->infile_num; infile_num_p++)
    {
        gzFile infile = gzopen((infile_stat->organized_infile_tab)[infile_num_p].fpath, "r");
        if (!infile)
            err(errno, "%s(): Cannot open file %s", __func__, (infile_stat->organized_infile_tab)[infile_num_p].fpath);
        kseq_t *seq = kseq_init(infile);

        while (kseq_read(seq) >= 0)
        {
            char *read = malloc(seq->seq.l + 1);
            if (!read)
                err(errno, "%dth genome malloc failed", infile_num_p);
            memcpy(read, seq->seq.s, seq->seq.l);
            read[seq->seq.l] = '\0';
            char **add = &read;
            vector_push(&all_reads, add);
            gseq_nums[(infile_num_p % batch_size) + 1]++; // gseq_nums[infile_num_p+1]++;
        }
        kseq_destroy(seq);
        gzclose(infile);
        if ((infile_num_p < infile_stat->infile_num - 1) && infile_num_p % batch_size < batch_size - 1)
            continue;

        int num = infile_num_p % batch_size + 1;
        for (int i = 0; i < num; i++)
            gseq_nums[i + 1] += gseq_nums[i];
#pragma omp parallel for num_threads(sketch_opt_val->p)
        for (uint32_t i = 0; i < num; i++)
        {
            khash_t(sort64) *h = kh_init(sort64);
            for (uint32_t j = gseq_nums[i]; j < gseq_nums[i + 1]; j++)
            {
                char *s = *(char **)vector_get(&all_reads, j);
                int len = strlen(s);
                if (len < klen)
                    continue;
                int base = 0;
                uint64_t tuple, crvstuple, unituple, basenum, unictx;
                for (int pos = 0; pos < len; pos++)
                {
                    if (Basemap[(unsigned short)s[pos]] == DEFAULT)
                    {
                        base = 0;
                        continue;
                    }
                    basenum = Basemap[(unsigned short)s[pos]];
                    tuple = ((tuple << 2) | basenum);
                    crvstuple = ((crvstuple >> 2) | ((basenum ^ 3LLU) << len_mv));
                    if (++base < klen)
                        continue;

                    unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
                    unictx = unituple & ctxmask;
                    // printf("%lx\t%lx\t%lx\t%lx\t%lx\t%lx\n", unituple, unictx, tuple,tuple & ctxmask, crvstuple, crvstuple & ctxmask);
                    if (SKETCH_HASH(unictx) > FILTER)
                        continue;

                    int ret;
                    khint_t key = kh_put(sort64, h, uint64kmer2generic_ctxobj(unituple & tupmask), &ret);
                    if (ret)
                        kh_value(h, key) = 1;
                } // nt pos loop
                free(s);
            } // seq j loop
            SortedKV_Arrays_t lco_ab = gpt_sort_khash_u64(h);
            // remove context with conflict object
            if (!sketch_opt_val->conflict)
                remove_ctx_with_conflict_obj(&lco_ab, Bitslen.obj);
            // may filter n for lco_ab
            batch_sketches[i] = lco_ab.keys;
            sketch_index[infile_num_p - num + i + 2] = lco_ab.len;
            kh_destroy(sort64, h);
            free(lco_ab.values);
        } // genome i loop

        for (uint32_t i = 0; i < num; i++)
        {
            fwrite(batch_sketches[i], sizeof(uint64_t), sketch_index[infile_num_p - num + i + 2], comb_sketch_fp);
            free(batch_sketches[i]);
        }

        memset(gseq_nums, 0, (batch_size + 1) * sizeof(int));
        vector_free(&all_reads);
        vector_init(&all_reads, sizeof(char *));

        printf("\r%d/%d genomes processed\n", infile_num_p + 1, infile_stat->infile_num);

    } // infile_num_p loop

    for (int i = 0; i < infile_stat->infile_num; i++)
        sketch_index[i + 1] += sketch_index[i];
    write_to_file(format_string("%s/%s", sketch_opt_val->outdir, idx_sketch_suffix), sketch_index, (infile_stat->infile_num + 1) * sizeof(sketch_index[0]));

    fclose(comb_sketch_fp);
    vector_free(&all_reads);
    free_all(sketch_index, gseq_nums, batch_sketches, NULL);
    write_sketch_stat(sketch_opt_val->outdir, infile_stat);
}

// for kssd3 ani use directly
simple_sketch_t *old_simple_genomes2mem2sortedctxobj64_mem(infile_tab_t *infile_stat, int drfold)
{ // only for coden pattern sketching*
    FILTER = UINT32_MAX >> drfold;
    uint32_t len_mv = 2 * klen - 2;

    uint64_t *sketch_index = calloc(infile_stat->infile_num + 1, sizeof(uint64_t));

    int *gseq_nums = calloc(infile_stat->infile_num + 1, sizeof(int));

    Vector all_reads;
    vector_init(&all_reads, sizeof(char *));

    int num = infile_stat->infile_num;
    uint64_t **batch_sketches = malloc(num * sizeof(uint64_t *));

    for (int infile_num_p = 0; infile_num_p < num; infile_num_p++)
    {
        gzFile infile = gzopen((infile_stat->organized_infile_tab)[infile_num_p].fpath, "r");
        if (!infile)
            err(errno, "%s(): Cannot open file %s", __func__, (infile_stat->organized_infile_tab)[infile_num_p].fpath);
        kseq_t *seq = kseq_init(infile);

        while (kseq_read(seq) >= 0)
        {
            char *read = malloc(seq->seq.l + 1);
            if (!read)
                err(errno, "%dth genome malloc failed", infile_num_p);
            memcpy(read, seq->seq.s, seq->seq.l);
            read[seq->seq.l] = '\0';
            char **add = &read;
            vector_push(&all_reads, add);
            gseq_nums[(infile_num_p) + 1]++; // gseq_nums[infile_num_p+1]++;
        }
        kseq_destroy(seq);
        gzclose(infile);

        if (infile_num_p < num - 1)
            continue;

        for (int i = 0; i < num; i++)
            gseq_nums[i + 1] += gseq_nums[i];

#pragma omp parallel for num_threads(2)
        for (uint32_t i = 0; i < num; i++)
        {

            khash_t(sort64) *h = kh_init(sort64);
            for (uint32_t j = gseq_nums[i]; j < gseq_nums[i + 1]; j++)
            {
                char *s = *(char **)vector_get(&all_reads, j);
                int len = strlen(s);
                if (len < klen)
                    continue;
                int base = 0;
                uint64_t tuple, crvstuple, unituple, basenum, unictx;
                for (int pos = 0; pos < len; pos++)
                {
                    if (Basemap[(unsigned short)s[pos]] == DEFAULT)
                    {
                        base = 0;
                        continue;
                    }
                    basenum = Basemap[(unsigned short)s[pos]];
                    tuple = ((tuple << 2) | basenum);
                    crvstuple = ((crvstuple >> 2) | ((basenum ^ 3LLU) << len_mv));
                    if (++base < klen)
                        continue;

                    unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
                    unictx = unituple & ctxmask;
                    if (SKETCH_HASH(unictx) > FILTER)
                        continue;

                    int ret;
                    khint_t key = kh_put(sort64, h, reorder_unituple_by_coden_pattern64(unituple & tupmask), &ret);
                    if (ret)
                        kh_value(h, key) = 1;
                } // nt pos loop
                free(s);
            } // seq j loop

            SortedKV_Arrays_t lco_ab = sort_khash_u64(h);

            // remove context with conflict object
            remove_ctx_with_conflict_obj(&lco_ab, Bitslen.obj);

            // may filter n for lco_ab
            batch_sketches[i] = lco_ab.keys;
            sketch_index[i + 1] = lco_ab.len;

            kh_destroy(sort64, h);
            free(lco_ab.values);
        } // genome i loop

    } // infile_num_p loop
    for (int i = 0; i < num; i++)
        sketch_index[i + 1] += sketch_index[i];

    uint64_t *combined_sketch = malloc(sizeof(uint64_t) * sketch_index[num]);
    if (!combined_sketch)
        err(EXIT_FAILURE, "combined_sketch allocation failed");

    for (int i = 0; i < num; i++)
    {
        memcpy(combined_sketch + sketch_index[i], batch_sketches[i], sizeof(uint64_t) * (sketch_index[i + 1] - sketch_index[i]));
        free(batch_sketches[i]);
    }

    vector_free(&all_reads);
    free_all(gseq_nums, batch_sketches, NULL);

    simple_sketch_t *return_sketch = malloc(sizeof(simple_sketch_t));
    if (!return_sketch)
        err(EXIT_FAILURE, "return_sketch allocation failed");

    return_sketch->comb_sketch = combined_sketch;
    return_sketch->sketch_index = sketch_index;
    return_sketch->infile_num = num;

    return return_sketch;
}

// drop-in replacement for read_genomes2mem2sortedctxobj64()
// keeps identical interface and output semantics
void test_read_genomes2mem2sortedctxobj64(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat, int batch_size)
{
    uint8_t nobjbits = Bitslen.obj;
    const uint32_t len_mv = (uint32_t)(2 * klen - 2); // shift for reverse rolling
    uint64_t *sketch_index = (uint64_t *)calloc((size_t)infile_stat->infile_num + 1, sizeof(uint64_t));
    if (!sketch_index)
        err(errno, "%s(): OOM sketch_index", __func__);

    // reuse across batches; only batch-local slots are touched
    uint64_t **batch_sketches = (uint64_t **)malloc((size_t)batch_size * sizeof(uint64_t *));
    if (!batch_sketches)
        err(errno, "%s(): OOM batch_sketches", __func__);

    // single combined output file, fully buffered
    FILE *comb_sketch_fp = fopen(format_string("%s/%s", sketch_opt_val->outdir, combined_sketch_suffix), "wb");
    if (!comb_sketch_fp)
        err(errno, "%s() open file error: %s/%s", __func__, sketch_opt_val->outdir, combined_sketch_suffix);
    // 8 MB buffered writes improve throughput on large outputs
    (void)setvbuf(comb_sketch_fp, NULL, _IOFBF, 8u << 20);
    //
    const int nfiles = infile_stat->infile_num;
    for (int batch_start = 0; batch_start < nfiles; batch_start += batch_size)
    {
        const int batch_end = (batch_start + batch_size <= nfiles) ? (batch_start + batch_size) : nfiles;
        const int this_batch = batch_end - batch_start;

        // lengths per genome in this batch
        uint64_t *lens = (uint64_t *)calloc((size_t)this_batch, sizeof(uint64_t));
        if (!lens)
            err(errno, "%s(): OOM lens", __func__);

// Process one genome per thread; no staging of reads in RAM.
#pragma omp parallel for num_threads(sketch_opt_val->p) schedule(dynamic, 1)
        for (int bi = 0; bi < this_batch; ++bi)
        {
            const int file_idx = batch_start + bi;
            const char *fpath = infile_stat->organized_infile_tab[file_idx].fpath;

            gzFile infile = gzopen(fpath, "r");
            if (!infile)
                err(errno, "%s(): Cannot open file %s", __func__, fpath);

            // enlarge zlib internal buffer to reduce syscalls
            (void)gzbuffer(infile, 1u << 20); // 1 MB

            kseq_t *seq = kseq_init(infile);
            if (!seq)
                err(errno, "%s(): kseq_init failed on %s", __func__, fpath);

            khash_t(sort64) *h = kh_init(sort64);
            if (!h)
                err(errno, "%s(): kh_init sort64 OOM", __func__);

            // Heuristic reserve: assume ~1/8 of bases produce valid kept k-mers after filtering.
            // Tweak if you know your sampling rate (c) to reduce rehash.
            kh_resize(sort64, h, 1u << 15); // start with 32K; grows as needed

            while (kseq_read(seq) >= 0)
            {
                const char *s = seq->seq.s;
                const int len = (int)seq->seq.l;
                if (len < klen)
                    continue;
                // rolling 2-bit canonical
                uint64_t tuple, crv; // forward and reverse-complement
                int base = 0;

                for (int pos = 0; pos < len; ++pos)
                {

                    const int bmap = Basemap[(unsigned char)s[pos]];
                    if (bmap == DEFAULT)
                    { // non-ACGT → reset window
                        base = 0;
                        continue;
                    }

                    const uint64_t b2 = (uint64_t)bmap;
                    tuple = (tuple << 2) | b2;
                    crv = (crv >> 2) | ((b2 ^ 3ull) << len_mv);
                    if (++base < klen)
                        continue;

                    // canonical by ctx, same as your original:
                    // compare (tuple & ctxmask) vs (crv & ctxmask)
                    const uint64_t t_ctx = (tuple & ctxmask);
                    const uint64_t r_ctx = (crv & ctxmask);
                    uint64_t unituple = (t_ctx < r_ctx) ? tuple : crv;
                    const uint64_t unictx = unituple & ctxmask;

                    // sketching decision based on context
                    if (SKETCH_HASH(unictx) > FILTER)
                        continue;
                    unituple &= tupmask;
                    // Rearrange to context-object (same length) and store encoded k-mer (not a hash)
                    // Your helper packs it; unchanged API:
                    // const uint64_t ctxobj = uint64kmer2generic_ctxobj(unituple & tupmask);
                    const uint64_t ctxobj = make_ctxobj(unituple, tupmask, ctxmask, nobjbits);

                    int ret;
                    khint_t key = kh_put(sort64, h, ctxobj, &ret);
                    if (ret)
                        kh_value(h, key) = 1; // presence; keep as in original
                } // end bases of this read
            } // end reads

            // finalize this genome: sort + (optionally) remove conflicts
            SortedKV_Arrays_t lco_ab = gpt_sort_khash_u64(h); // sort_khash_u64(h);
            if (!sketch_opt_val->conflict)
                remove_ctx_with_conflict_obj(&lco_ab, nobjbits);

            batch_sketches[bi] = lco_ab.keys; // will be written by main thread
            lens[bi] = lco_ab.len;

            kh_destroy(sort64, h);
            free(lco_ab.values);

            kseq_destroy(seq);
            gzclose(infile);
        } // end omp per-genome

        // Write this batch in file order, build sketch_index
        for (int bi = 0; bi < this_batch; ++bi)
        {
            const int global_i = batch_start + bi;
            const uint64_t n = lens[bi];

            // index prefix sums are built after all batches (like your original),
            // here we record per-genome lengths at sketch_index[i+1]
            sketch_index[global_i + 1] = n;

            if (n)
            {
                fwrite(batch_sketches[bi], sizeof(uint64_t), n, comb_sketch_fp);
                free(batch_sketches[bi]);
                batch_sketches[bi] = NULL;
            }
        }

        free(lens);

        printf("\r%d/%d genomes processed\n", batch_end, nfiles);
        fflush(stdout);
    } // end batches

    // prefix sum → absolute offsets (same as your original)
    for (int i = 0; i < nfiles; ++i)
        sketch_index[i + 1] += sketch_index[i];

    write_to_file(format_string("%s/%s", sketch_opt_val->outdir, idx_sketch_suffix),
                  sketch_index,
                  (size_t)(nfiles + 1) * sizeof(sketch_index[0]));

    fclose(comb_sketch_fp);
    free(batch_sketches);
    free(sketch_index);

    write_sketch_stat(sketch_opt_val->outdir, infile_stat);
}


// === unified helpers + both modes (A and B) ===
// ---- conflict filters provided by your code:
// void remove_ctx_with_conflict_obj(SortedKV_Arrays_t *kv, uint32_t n_obj_bits);
// void remove_ctx_with_conflict_obj_noabund(uint64_t *keys, size_t *len_io, uint32_t n_obj_bits);
#include <omp.h>
typedef struct{char *s; int l;} read_span_t;
// ---- In-file parallel helpers (used only when #files < threads) ----

// Hot loop: push kept ctxobj into a vector (no hashing)
static inline void sketch_read_into_vec(const char *restrict s, int len, u64vec *restrict vec,
                                        uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits, uint32_t klen)
{
    if (len < (int)klen) return;
    const uint32_t len_mv = (uint32_t)(2 * klen - 2);

    uint64_t tuple = 0, crv = 0;
    int base = 0;

    for (int pos = 0; pos < len; ++pos)
    {
        const int bmap = Basemap[(unsigned char)s[pos]];
        if (unlikely(bmap == DEFAULT))
        {
            base = 0;
            tuple = 0;
            crv = 0;
            continue;
        }

        const uint64_t b2 = (uint64_t)bmap;
        tuple = (tuple << 2) | b2;
        crv = (crv >> 2) | ((b2 ^ 3ull) << len_mv);
        if (unlikely(++base < (int)klen)) continue;

        const uint64_t t_ctx = (tuple & ctxmask);
        const uint64_t r_ctx = (crv & ctxmask);
        const int use_fwd = (t_ctx < r_ctx);
        const uint64_t unictx = use_fwd ? t_ctx : r_ctx;

        if (unlikely(SKETCH_HASH(unictx) > FILTER)) continue;

        const uint64_t unituple = (use_fwd ? tuple : crv) & tupmask;
        const uint64_t ctxobj = make_ctxobj(unituple, tupmask, ctxmask, nobjbits);
        v_push(vec, ctxobj);
    }
}

static inline void shrink_thread_vec(u64vec *vec, uint32_t **counts_out)
{
    // vec->a contains keys; produce counts_out aligned to deduped keys
    if (vec->n == 0) { *counts_out = NULL; return;}
    if (vec->n == 1){ // fast path
        *counts_out = (uint32_t *)malloc(sizeof(uint32_t));
        if (*counts_out) (*counts_out)[0] = 1;
        return;
    }
    radix_sort_u64(vec->a, vec->n);
    vec->n = dedup_with_counts(vec->a, vec->n, counts_out);
    if (!*counts_out)
    { // fallback if alloc failed inside helper
        *counts_out = (uint32_t *)malloc(vec->n * sizeof(uint32_t));
        if (*counts_out)
            for (size_t i = 0; i < vec->n; ++i)
                (*counts_out)[i] = 1u;
    }
}


static inline SortedKV_Arrays_t build_kv_from_vec(u64vec *vec, bool has_abundance)
{
    SortedKV_Arrays_t kv = (SortedKV_Arrays_t){0};
    if (vec->n == 0) return kv;

    radix_sort_u64(vec->a, vec->n);

    if (has_abundance) {
        uint32_t *counts = NULL;
        vec->n = dedup_with_counts(vec->a, vec->n, &counts);
        if (!counts) {
            counts = (uint32_t*)malloc(vec->n * sizeof(uint32_t));
            if (!counts) return kv; // OOM; kv stays empty, caller handles
            for (size_t i=0;i<vec->n;++i) counts[i] = 1u;
        }
        kv.keys   = (uint64_t*)malloc(vec->n * sizeof(uint64_t));
        kv.values = (uint32_t*)malloc(vec->n * sizeof(uint32_t));
        if (!kv.keys || !kv.values) { free(kv.keys); free(kv.values); free(counts); return (SortedKV_Arrays_t){0}; }
        memcpy(kv.keys, vec->a, vec->n * sizeof(uint64_t));
        memcpy(kv.values, counts, vec->n * sizeof(uint32_t));
        kv.len = vec->n;
        free(counts);
    } else {
        vec->n = dedup_sorted_uint64(vec->a, vec->n);
        kv.keys = (uint64_t*)malloc(vec->n * sizeof(uint64_t));
        if (!kv.keys) return (SortedKV_Arrays_t){0};
        memcpy(kv.keys, vec->a, vec->n * sizeof(uint64_t));
        kv.values = NULL;
        kv.len    = vec->n;
    }
    return kv;
}

// ---------- unified FASTQ loader (gz or plain) → u64vec ----------
// Requires:
//   - sketch_read_into_vec(const char*, int, u64vec*, uint64_t, uint64_t, uint8_t, uint32_t)
//   - globals/params: ctxmask, tupmask, Bitslen.obj, klen
//   - zlib/kseq headers linked (for gz path)
// Build: add -fopenmp if you parallelize elsewhere; this function itself is single-threaded.

// small helper: check suffix (case-sensitive)
static inline int has_suffix(const char *s, const char *suf) {
    size_t ns=0, nf=0; while (s[ns]) ++ns; while (suf[nf]) ++nf;
    return (nf<=ns) && (memcmp(s+ns-nf, suf, nf)==0);
}

// =====================================================
// Mode-A helpers: single-vector collectors (gz / plain)
// =====================================================

// Single-vector collector for gz FASTQ
static inline void load_fastx_into_single_vec(char *path, u64vec *vec, uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits,uint32_t klen)
{
    gzFile in = gzopen(path, "r");
    if (!in) err(errno, "%s(): Cannot open %s", __func__, path);
    (void)gzbuffer(in, 4u<<20);

    kseq_t *seq = kseq_init(in);
    if (!seq) err(errno, "%s(): kseq_init %s", __func__, path);

    while (kseq_read(seq) >= 0) {
        if (seq->seq.l) {
            sketch_read_into_vec(seq->seq.s, (int)seq->seq.l, vec,
                                 ctxmask, tupmask, nobjbits, klen);
        }
    }
    kseq_destroy(seq);
    gzclose(in);
}

// Single-vector collector for plain FASTQ via mmap (no copies)
static inline void load_fastq_plain_mmap_into_single_vec(char *path, u64vec *vec, uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits, uint32_t klen)
{
    int fd = open(path, O_RDONLY);
    if (fd < 0) err(errno, "%s(): open %s", __func__, path);
#ifdef POSIX_FADV_SEQUENTIAL
    posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    struct stat st;
    if (fstat(fd, &st) != 0) err(errno, "%s(): fstat %s", __func__, path);
    const size_t fsz = (size_t)st.st_size;
    if (fsz == 0) { close(fd); return; }

    char *base = (char*)mmap(NULL, fsz, PROT_READ, MAP_PRIVATE, fd, 0);
    if (base == MAP_FAILED) err(errno, "%s(): mmap %s", __func__, path);
#ifdef MADV_SEQUENTIAL
    madvise(base, fsz, MADV_SEQUENTIAL);
#endif
#ifdef MADV_WILLNEED
    madvise(base, fsz, MADV_WILLNEED);
#endif

    // Walk lines; every 2nd of each 4 is the sequence
    size_t line_start = 0;
    unsigned line_mod = 0; // 0=@hdr,1=seq,2=+,3=qual
    for (size_t i = 0; i < fsz; ++i) {
        if (base[i] != '\n') continue;
        size_t L = i - line_start;
        if (L && base[i-1] == '\r') --L; // trim CR
        if (line_mod == 1 && L) {
            sketch_read_into_vec(base + line_start, (int)L, vec,
                                 ctxmask, tupmask, nobjbits, klen);
        }
        line_start = i + 1;
        line_mod = (line_mod + 1) & 3;
    }
    // (Optional) handle last unterminated line — FASTQ typically ends with '\n'

    munmap(base, fsz);
    close(fd);
}

// Tiny wrapper to choose gz vs plain
static inline void load_genome_into_single_vec(char *path,u64vec *vec,uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits,uint32_t klen)
{
    if (isCompressfile(path) || isOK_fmt_infile(path, fasta_fmt, FAS_FMT_SZ))
        load_fastx_into_single_vec(path, vec, ctxmask, tupmask, nobjbits, klen);
    else if (isOK_fmt_infile(path, fastq_fmt, FQ_FMT_SZ))
        load_fastq_plain_mmap_into_single_vec(path, vec, ctxmask, tupmask, nobjbits, klen);
    else  err(errno, "%s(): %s is not accept format(.fasta,.fastq)",__func__, path);    
}

// =====================================================
// ---- Mode A: per-genome parallel (vector + counts) --
// =====================================================
static void sketch_many_files_in_parallel(sketch_opt_t *opt, infile_tab_t *tab, int batch_size)
{
    const int nfiles = tab->infile_num;
    const int kmerocrs = opt->kmerocrs;              
    const bool has_abundance = opt->abundance ;

    FILE *comb = fopen(format_string("%s/%s", opt->outdir, combined_sketch_suffix), "wb");
    if (!comb) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_sketch_suffix);
    setvbuf(comb, NULL, _IOFBF, 8u << 20);

    FILE *comb_ab = NULL;
    if (has_abundance) {
        comb_ab = fopen(format_string("%s/%s", opt->outdir, combined_ab_suffix), "wb");
        if (!comb_ab) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_ab_suffix);
        setvbuf(comb_ab, NULL, _IOFBF, 8u << 20);
    }

    uint64_t *sketch_index = (uint64_t *)calloc((size_t)nfiles + 1, sizeof(uint64_t));
    if (!sketch_index) err(errno, "%s(): OOM index", __func__);

    // Process in batches of files to cap memory if needed
    for (int batch_start = 0; batch_start < nfiles; batch_start += batch_size) {
        const int batch_end  = (batch_start + batch_size <= nfiles) ? (batch_start + batch_size) : nfiles;
        const int this_batch = batch_end - batch_start;

        uint64_t *lens = (uint64_t *)calloc((size_t)this_batch, sizeof(uint64_t));
        if (!lens) err(errno, "%s(): OOM lens", __func__);

        uint64_t **batch_keys = (uint64_t **)calloc((size_t)this_batch, sizeof(uint64_t *));
        if (!batch_keys) err(errno, "%s(): OOM batch_keys", __func__);

        uint32_t **batch_vals = NULL;
        if (has_abundance) {
            batch_vals = (uint32_t **)calloc((size_t)this_batch, sizeof(uint32_t *));
            if (!batch_vals) err(errno, "%s(): OOM batch_vals", __func__);
        }

        // One thread per file (parallel across files)
        #pragma omp parallel for num_threads(opt->p) schedule(dynamic, 1)
        for (int bi = 0; bi < this_batch; ++bi) {
            const int file_idx = batch_start + bi;
            char *fpath  = tab->organized_infile_tab[file_idx].fpath;

            u64vec vec; v_init(&vec, 1u << 15);
            load_genome_into_single_vec(fpath, &vec, ctxmask, tupmask, Bitslen.obj, klen);

            if (vec.n) {
                SortedKV_Arrays_t kv = build_kv_from_vec(&vec, has_abundance || (kmerocrs > 1) );
                if(kmerocrs > 1) // filter by min ocrs
                    filter_n_SortedKV_Arrays(&kv, kmerocrs);
                if (!opt->conflict){
                    if(has_abundance) remove_ctx_with_conflict_obj(&kv, Bitslen.obj);
                    else remove_ctx_with_conflict_obj_noabund (kv.keys, &(kv.len), Bitslen.obj);
                }
                batch_keys[bi] = kv.keys;
                lens[bi]       = kv.len;
                if (has_abundance) batch_vals[bi] = kv.values;
                else               free(kv.values);
            } else {
                lens[bi] = 0;
            }
            v_free(&vec);
        }

        // write batch in order
        for (int bi = 0; bi < this_batch; ++bi) {
            const int gi = batch_start + bi;
            sketch_index[gi + 1] = lens[bi];
            if (lens[bi]) {
                fwrite(batch_keys[bi], sizeof(uint64_t), lens[bi], comb);
                free(batch_keys[bi]);
                if (has_abundance) {
                    fwrite(batch_vals[bi], sizeof(uint32_t), lens[bi], comb_ab);
                    free(batch_vals[bi]);
                }
            }
        }

        free(batch_keys);
        free(lens);
        if (has_abundance) free(batch_vals);

        fprintf(stderr, "\r%d/%d genomes processed", batch_end, nfiles);
        fflush(stderr);
    }

    // prefix-sum index and finish
    for (int i = 0; i < nfiles; ++i) sketch_index[i + 1] += sketch_index[i];
    write_to_file(format_string("%s/%s", opt->outdir, idx_sketch_suffix),
                  sketch_index, (size_t)(nfiles + 1) * sizeof(uint64_t));
    fclose(comb);
    if (has_abundance) fclose(comb_ab);
    free(sketch_index);
    write_sketch_stat(opt->outdir, tab);
}

// ============================
// Unified helpers + Mode-B
// ============================


// C helper for finding '\n' between [off, lim)
static inline char *find_nl(const char *base, size_t off, size_t lim) {
    if (off >= lim) return NULL;
    return (char*)memchr(base + off, '\n', lim - off);
}

// ---------- shrink helpers ----------
static inline unsigned clamp_bucket_bits_from_total(size_t total_distinct){
    // target ~4k entries per bucket; clamp to [8..14]
    size_t target_buckets = (total_distinct / 4096) ? (total_distinct / 4096) : 256;
    unsigned BITS = 8;
    while (((size_t)1 << BITS) < target_buckets && BITS < 14) ++BITS;
    return BITS;
}

static inline size_t shrink_kv_inplace_u64_u32(uint64_t *k, uint32_t *v, size_t n){
    if (!n) return 0;
    size_t w=0;
    for (size_t r=1; r<n; ++r){
        if (k[r]==k[w]) v[w] += v[r];
        else { ++w; k[w]=k[r]; v[w]=v[r]; }
    }
    return w+1;
}

// ---------- collectors (no finalization here) ----------
static void collect_gz_into_vectors(const char *path, int nth, int nthreads_for_omp,
                                    u64vec *V_thr,
                                    uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits, uint32_t klen,
                                    int BATCH_READS)
{
    gzFile in = gzopen(path, "r");
    if (!in) err(errno, "%s(): Cannot open %s", __func__, path);
    (void)gzbuffer(in, 4u<<20);
    kseq_t *seq = kseq_init(in);
    if (!seq) err(errno, "%s(): kseq_init %s", __func__, path);

    read_span_t *R = (read_span_t*)malloc((size_t)BATCH_READS*sizeof(*R));
    if (!R) err(errno, "%s(): OOM batch reads", __func__);
    int rcnt=0;

    int ret;
    do{
        ret = kseq_read(seq);
        if (ret >= 0){
            char *buf = (char*)malloc(seq->seq.l + 1);
            if (!buf) err(errno, "%s(): OOM read buf", __func__);
            memcpy(buf, seq->seq.s, seq->seq.l);
            buf[seq->seq.l] = '\0';
            R[rcnt].s = buf; R[rcnt].l = (int)seq->seq.l;
            ++rcnt;
        }
        if (rcnt >= BATCH_READS || ret < 0){
            #pragma omp parallel for if (nth>1) num_threads(nthreads_for_omp) schedule(static)
            for (int ri=0; ri<rcnt; ++ri){
                int tid = 0;
                #ifdef _OPENMP
                tid = omp_get_thread_num();
                #endif
                sketch_read_into_vec(R[ri].s, R[ri].l, &V_thr[tid],
                                     ctxmask, tupmask, nobjbits, klen);
            }
            for (int i=0;i<rcnt;++i) free(R[i].s);
            rcnt=0;
        }
    } while (ret >= 0);
    free(R);
    kseq_destroy(seq);
    gzclose(in);
}

static void collect_plain_mmap_into_vectors(const char *path, int nth, int nthreads_for_omp,
                                            u64vec *V_thr,
                                            uint64_t ctxmask, uint64_t tupmask, uint8_t nobjbits, uint32_t klen)
{
    int fd = open(path, O_RDONLY);
    if (fd < 0) err(errno, "%s(): open %s", __func__, path);
#ifdef POSIX_FADV_SEQUENTIAL
    posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    struct stat st;
    if (fstat(fd, &st) != 0) err(errno, "%s(): fstat %s", __func__, path);
    const size_t fsz = (size_t)st.st_size;
    if (fsz == 0) { close(fd); return; }

    char *base = (char*)mmap(NULL, fsz, PROT_READ, MAP_PRIVATE, fd, 0);
    if (base == MAP_FAILED) err(errno, "%s(): mmap %s", __func__, path);
#ifdef MADV_SEQUENTIAL
    madvise(base, fsz, MADV_SEQUENTIAL);
#endif
#ifdef MADV_WILLNEED
    madvise(base, fsz, MADV_WILLNEED);
#endif

    const size_t chunk = (fsz + (size_t)nth - 1) / (size_t)nth;

    #pragma omp parallel for if (nth>1) num_threads(nthreads_for_omp) schedule(static)
    for (int tid=0; tid<nth; ++tid){
        size_t start = (size_t)tid * chunk;
        size_t end   = (tid == nth-1) ? fsz : start + chunk;
        if (start >= fsz) { V_thr[tid].n = 0; continue; }
        if (start > 0) { while (start < end && base[start-1] != '\n') ++start; }

        size_t cursor = start;
        while (cursor < end) {
            if (base[cursor] == '@') {
                char *nl1 = find_nl(base, cursor, fsz);                 if (!nl1) break;
                size_t seq_beg = (size_t)(nl1 - base) + 1;
                char *nl2 = find_nl(base, seq_beg, fsz);                if (!nl2) break;
                size_t plus_beg = (size_t)(nl2 - base) + 1;
                char *nl3 = find_nl(base, plus_beg, fsz);               if (!nl3) break;
                size_t qual_beg = (size_t)(nl3 - base) + 1;
                char *nl4 = find_nl(base, qual_beg, fsz);               if (!nl4) break;

                if (cursor >= start) {
                    size_t L = (size_t)(nl2 - (base + seq_beg));
                    if (L && base[seq_beg + L - 1] == '\r') --L;
                    if (L) {
                        sketch_read_into_vec(base + seq_beg, (int)L, &V_thr[tid],
                                             ctxmask, tupmask, nobjbits, klen);
                    }
                }
                cursor = (size_t)(nl4 - base) + 1;
                if (cursor >= end && tid != nth-1) break;
            } else {
                char *nl = find_nl(base, cursor, end);
                if (!nl) break;
                cursor = (size_t)(nl - base) + 1;
            }
        }
    }

    munmap(base, fsz);
    close(fd);
}

// ---------- shared finalizer (shrink + bucketed merge) ----------
static SortedKV_Arrays_t finalize_vectors_bucketed(u64vec *V_thr, int nth, int nthreads_for_omp)
{
    SortedKV_Arrays_t out = (SortedKV_Arrays_t){0};
    uint32_t **C_thr = (uint32_t**)calloc((size_t)nth, sizeof(uint32_t*));
    if (!C_thr) return out;

    size_t pre_elems=0; for (int t=0;t<nth;++t) pre_elems += V_thr[t].n;

    #pragma omp parallel for if (nth>1 && pre_elems >= (1u<<20)) num_threads(nthreads_for_omp) schedule(dynamic)
    for (int t=0; t<nth; ++t) shrink_thread_vec(&V_thr[t], &C_thr[t]);

    size_t total_distinct=0; for (int t=0;t<nth;++t) total_distinct += V_thr[t].n;
    if (total_distinct == 0){
        for (int t=0;t<nth;++t){ v_free(&V_thr[t]); free(C_thr[t]); }
        free(V_thr); free(C_thr);
        return out;
    }

    const unsigned BITS = clamp_bucket_bits_from_total(total_distinct);
    const size_t   NB   = (size_t)1u << BITS;

    size_t *bucket_totals = (size_t*)calloc(NB, sizeof(size_t));
    size_t **thr_counts   = (size_t**)malloc((size_t)nth * sizeof(size_t*));
    if (!bucket_totals || !thr_counts) { free(bucket_totals); free(thr_counts); goto CLEAN_EMPTY; }

    for (int t=0;t<nth;++t){
        thr_counts[t] = (size_t*)calloc(NB, sizeof(size_t));
        if (!thr_counts[t]) { for (int j=0;j<t;++j) free(thr_counts[j]); free(thr_counts); free(bucket_totals); goto CLEAN_EMPTY; }
    }

    #pragma omp parallel for num_threads(nthreads_for_omp) schedule(static)
    for (int t=0; t<nth; ++t){
        uint64_t *A = V_thr[t].a; size_t n = V_thr[t].n;
        for (size_t i=0;i<n;++i){
            unsigned b = (unsigned)(A[i] >> (64 - BITS));
            ++thr_counts[t][b];
        }
    }
    for (size_t b=0;b<NB;++b){
        size_t sum=0; for (int t=0;t<nth;++t) sum += thr_counts[t][b];
        bucket_totals[b] = sum;
    }

    size_t *bucket_off = (size_t*)malloc((NB+1)*sizeof(size_t));
    if (!bucket_off) { for (int t=0;t<nth;++t) free(thr_counts[t]); free(thr_counts); free(bucket_totals); goto CLEAN_EMPTY; }
    bucket_off[0]=0; for (size_t b=0;b<NB;++b) bucket_off[b+1] = bucket_off[b] + bucket_totals[b];
    const size_t TOTAL = bucket_off[NB];

    uint64_t *B_keys = (uint64_t*)malloc(TOTAL*sizeof(uint64_t));
    uint32_t *B_vals = (uint32_t*)malloc(TOTAL*sizeof(uint32_t));
    if (!B_keys || !B_vals) { free(B_keys); free(B_vals); free(bucket_off); for (int t=0;t<nth;++t) free(thr_counts[t]); free(thr_counts); free(bucket_totals); goto CLEAN_EMPTY; }

    size_t **thr_write = (size_t**)malloc((size_t)nth * sizeof(size_t*));
    if (!thr_write) { free(B_keys); free(B_vals); free(bucket_off); for (int t=0;t<nth;++t) free(thr_counts[t]); free(thr_counts); free(bucket_totals); goto CLEAN_EMPTY; }
    for (int t=0;t<nth;++t){
        thr_write[t] = (size_t*)malloc(NB*sizeof(size_t));
        if (!thr_write[t]) { for (int j=0;j<t;++j) free(thr_write[j]); free(thr_write); free(B_keys); free(B_vals); free(bucket_off); for (int j2=0;j2<nth;++j2) free(thr_counts[j2]); free(thr_counts); free(bucket_totals); goto CLEAN_EMPTY; }
    }
    for (size_t b=0;b<NB;++b){
        size_t off = bucket_off[b];
        for (int t=0;t<nth;++t){ thr_write[t][b] = off; off += thr_counts[t][b]; }
    }

    #pragma omp parallel for num_threads(nthreads_for_omp) schedule(static)
    for (int t=0; t<nth; ++t){
        uint64_t *Ak = V_thr[t].a; uint32_t *Av = C_thr[t]; size_t n = V_thr[t].n;
        size_t *pos = thr_write[t];
        for (size_t i=0;i<n;++i){
            uint64_t k = Ak[i];
            unsigned b = (unsigned)(k >> (64 - BITS));
            size_t p = pos[b]++;
            B_keys[p] = k;
            B_vals[p] = Av ? Av[i] : 1u;
        }
    }

    for (int t=0;t<nth;++t){ v_free(&V_thr[t]); free(C_thr[t]); free(thr_counts[t]); free(thr_write[t]); }
    free(V_thr); free(C_thr); free(thr_counts); free(thr_write);

    uint64_t **bk_keys = (uint64_t**)malloc(NB*sizeof(uint64_t*));
    uint32_t **bk_vals = (uint32_t**)malloc(NB*sizeof(uint32_t*));
    size_t    *bk_len  = (size_t*)   malloc(NB*sizeof(size_t));
    if (!bk_keys || !bk_vals || !bk_len) { free(bk_keys); free(bk_vals); free(bk_len); free(B_keys); free(B_vals); free(bucket_off); free(bucket_totals); goto CLEAN_EMPTY; }

    #pragma omp parallel for num_threads(nthreads_for_omp) schedule(dynamic)
    for (size_t b=0; b<NB; ++b){
        const size_t begin = bucket_off[b], end = bucket_off[b+1];
        const size_t n = end - begin;
        if (!n){ bk_keys[b]=NULL; bk_vals[b]=NULL; bk_len[b]=0; continue; }
        uint64_t *kseg = B_keys + begin; uint32_t *vseg = B_vals + begin;
        radix_sort_kv_u64(kseg, vseg, n);
        size_t m = shrink_kv_inplace_u64_u32(kseg, vseg, n);
        uint64_t *ko = (uint64_t*)malloc(m*sizeof(uint64_t));
        uint32_t *vo = (uint32_t*)malloc(m*sizeof(uint32_t));
        if (!ko || !vo) err(EXIT_FAILURE, "%s(): OOM bucket kv out", __func__);
        memcpy(ko, kseg, m*sizeof(uint64_t));
        memcpy(vo, vseg, m*sizeof(uint32_t));
        bk_keys[b]=ko; bk_vals[b]=vo; bk_len[b]=m;
    }

    free(B_keys); free(B_vals); free(bucket_totals);

    size_t *out_off = (size_t*)malloc((NB+1)*sizeof(size_t));
    if (!out_off) { free(bk_keys); free(bk_vals); free(bk_len); free(bucket_off); goto CLEAN_EMPTY; }
    out_off[0]=0; for (size_t b=0;b<NB;++b) out_off[b+1] = out_off[b] + bk_len[b];
    const size_t M = out_off[NB];

    out.keys   = (uint64_t*)malloc(M*sizeof(uint64_t));
    out.values = (uint32_t*)malloc(M*sizeof(uint32_t));
    if (!out.keys || !out.values) { free(out.keys); free(out.values); out = (SortedKV_Arrays_t){0}; }

    #pragma omp parallel for num_threads(nthreads_for_omp) schedule(static)
    for (size_t b=0; b<NB; ++b){
        if (!bk_len[b]) continue;
        memcpy(out.keys  + out_off[b], bk_keys[b], bk_len[b]*sizeof(uint64_t));
        memcpy(out.values+ out_off[b], bk_vals[b], bk_len[b]*sizeof(uint32_t));
        free(bk_keys[b]); free(bk_vals[b]);
    }
    free(bk_keys); free(bk_vals); free(bk_len); free(bucket_off); free(out_off);

    out.len = M;
    return out;

CLEAN_EMPTY:
    for (int t=0;t<nth;++t){ v_free(&V_thr[t]); /* C_thr[t] may be NULL */ free(C_thr[t]); }
    free(V_thr); free(C_thr);
    return (SortedKV_Arrays_t){0};
}


static inline void write_index_payload_and_free(FILE *comb, FILE *comb_ab,
                bool write_ab, uint64_t *sketch_index_slot, SortedKV_Arrays_t *kv)
{
    *sketch_index_slot = kv->len;
    if (kv->len) {
        fwrite(kv->keys,   sizeof(uint64_t), kv->len, comb);
        if (write_ab) fwrite(kv->values, sizeof(uint32_t), kv->len, comb_ab);
    }
    free(kv->keys);   kv->keys   = NULL;
    free(kv->values); kv->values = NULL;
    kv->len = 0;
}

// =================== Mode B (in-file parallel; gz + mmap) ===================
static void sketch_few_files_with_intrafile_parallel(sketch_opt_t *opt, infile_tab_t *tab, int BATCH_READS)
{
    const int nfiles = tab->infile_num;
    const int kmerocrs = opt->kmerocrs;
    FILE *comb = fopen(format_string("%s/%s", opt->outdir, combined_sketch_suffix), "wb");
    if (!comb) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_sketch_suffix);
    setvbuf(comb, NULL, _IOFBF, 8u << 20);

    FILE *comb_ab = NULL;
    if (opt->abundance) {
        comb_ab = fopen(format_string("%s/%s", opt->outdir, combined_ab_suffix), "wb");
        if (!comb_ab) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_ab_suffix);
        setvbuf(comb_ab, NULL, _IOFBF, 8u<<20);
    }

    uint64_t *sketch_index = (uint64_t*)calloc((size_t)nfiles + 1, sizeof(uint64_t));
    if (!sketch_index) err(errno, "%s(): OOM sketch_index", __func__);

    for (int f=0; f<nfiles; ++f){
        const char *path = tab->organized_infile_tab[f].fpath;
        const int nth = (opt->p > 0 ? opt->p : 1);

        // per-thread accumulators
        u64vec *V_thr = (u64vec*)malloc((size_t)nth*sizeof(u64vec));
        if (!V_thr) err(errno, "%s(): OOM V_thr", __func__);
        for (int t=0;t<nth;++t) v_init(&V_thr[t], 1u<<15);

        // 1) collect (mmap for plain; kseq for .gz)
        if (isCompressfile((char*)path) || isOK_fmt_infile((char*)path, fasta_fmt, FAS_FMT_SZ))
            collect_gz_into_vectors(path, nth, opt->p, V_thr, ctxmask, tupmask, Bitslen.obj, klen, BATCH_READS);
        else // collect_plain_mmap_into_vectors() only for plain FASTQ not fast
            collect_plain_mmap_into_vectors(path, nth, opt->p, V_thr, ctxmask, tupmask, Bitslen.obj, klen);
        

        // 2) finalize once (shared shrink + bucketed merge)
        SortedKV_Arrays_t kv = finalize_vectors_bucketed(V_thr, nth, opt->p);
        if(kmerocrs > 1) // filter by min ocrs
                    filter_n_SortedKV_Arrays(&kv, kmerocrs);
        // 3) conflict filter + 4) write
        if (kv.len){
            if (!opt->conflict){
                    if(opt->abundance) remove_ctx_with_conflict_obj(&kv, Bitslen.obj);
                    else remove_ctx_with_conflict_obj_noabund (kv.keys, &(kv.len), Bitslen.obj);
            }
            write_index_payload_and_free(comb, (opt->abundance?comb_ab:NULL),
                                         opt->abundance, &sketch_index[f+1], &kv);
        } else {
            sketch_index[f+1] = 0;
        }

        fprintf(stderr, "\r%d/%d genomes processed", f+1, nfiles);
        fflush(stderr);
    }

    for (int i=0;i<nfiles;++i) sketch_index[i+1] += sketch_index[i];
    write_to_file(format_string("%s/%s", opt->outdir, idx_sketch_suffix),
                  sketch_index, (size_t)(nfiles + 1) * sizeof(uint64_t));
    fclose(comb);
    if (opt->abundance) fclose(comb_ab);
    free(sketch_index);
    write_sketch_stat(opt->outdir, tab);
}

// simple API for kssd3 ani use directly
// Faster, memory-only, keys-only sketcher (vector + radix + dedup + conflict filter)
// API kept identical to your original.

simple_sketch_t *simple_genomes2mem2sortedctxobj64_mem(infile_tab_t *infile_stat, int drfold)
{
    // ---- 1) configure sketch filter (same semantics as your original) ----
    FILTER = UINT32_MAX >> drfold;

    const int nfiles = infile_stat->infile_num;
    if (nfiles <= 0) {
        simple_sketch_t *ret = (simple_sketch_t*)calloc(1, sizeof(simple_sketch_t));
        if (!ret) err(EXIT_FAILURE, "%s(): OOM ret", __func__);
        ret->comb_sketch = NULL;
        ret->sketch_index = (uint64_t*)calloc(1, sizeof(uint64_t)); // [0] = 0
        ret->infile_num = 0;
        return ret;
    }

    // ---- 2) per-file outputs (keys only) ----
    uint64_t **per_keys = (uint64_t**)calloc((size_t)nfiles, sizeof(uint64_t*));
    uint64_t  *per_len  = (uint64_t*) calloc((size_t)nfiles, sizeof(uint64_t));
    if (!per_keys || !per_len) err(EXIT_FAILURE, "%s(): OOM per-file arrays", __func__);

    // ---- 3) parallel across files: load → vectorize (ctxobj) → sort → dedup → conflict-filter ----
    // Set number of threads via OMP env or your build flags (-fopenmp)
    #pragma omp parallel for schedule(dynamic,1) num_threads(nfiles)
    for (int f = 0; f < nfiles; ++f) {
        char *path = infile_stat->organized_infile_tab[f].fpath;

        // Collect kept ctxobj keys into a single vector for this file
        u64vec vec; v_init(&vec, 1u << 15); // heuristic starting capacity
        load_genome_into_single_vec(path, &vec, ctxmask, tupmask, Bitslen.obj, klen);
        if (vec.n == 0) {
            per_keys[f] = NULL;
            per_len[f]  = 0;
            v_free(&vec);
            continue;
        }

        // Sort + dedup (keys only)
        radix_sort_u64(vec.a, vec.n);
        vec.n = dedup_sorted_uint64(vec.a, vec.n);
        // Conflict filter (no abundance version) — keeps array sorted, shrinks length in-place
        remove_ctx_with_conflict_obj_noabund(vec.a, &vec.n, Bitslen.obj);

        // Transfer ownership: shrink to exact size and hand out
        uint64_t *outk = (uint64_t*)malloc(vec.n * sizeof(uint64_t));
        if (!outk) err(EXIT_FAILURE, "%s(): OOM outk", __func__);
        memcpy(outk, vec.a, vec.n * sizeof(uint64_t));
        per_keys[f] = outk;
        per_len[f]  = (uint64_t)vec.n;

        v_free(&vec);
    }

    // ---- 4) prefix-sum index & concatenate all into one combined sketch buffer ----
    uint64_t *sketch_index = (uint64_t*)calloc((size_t)nfiles + 1, sizeof(uint64_t));
    if (!sketch_index) err(EXIT_FAILURE, "%s(): OOM sketch_index", __func__);

    for (int i = 0; i < nfiles; ++i) sketch_index[i + 1] = sketch_index[i] + per_len[i];
    const uint64_t total_keys = sketch_index[nfiles];

    uint64_t *combined = NULL;
    if (total_keys) {
        combined = (uint64_t*)malloc((size_t)total_keys * sizeof(uint64_t));
        if (!combined) err(EXIT_FAILURE, "%s(): OOM combined", __func__);

        // copy each file's keys to its slot
        for (int i = 0; i < nfiles; ++i) {
            const uint64_t n = per_len[i];
            if (!n) continue;
            memcpy(combined + sketch_index[i], per_keys[i], (size_t)n * sizeof(uint64_t));
            free(per_keys[i]); // done with this piece
        }
    }
    free(per_keys);
    free(per_len);

    // ---- 5) build return object ----
    simple_sketch_t *ret = (simple_sketch_t*)malloc(sizeof(simple_sketch_t));
    if (!ret) err(EXIT_FAILURE, "%s(): OOM ret", __func__);
    ret->comb_sketch = combined;      // length = total_keys (may be 0)
    ret->sketch_index = sketch_index; // length = nfiles+1
    ret->infile_num = nfiles;
    return ret;
}



typedef struct {
    char *name;
    char *seq;
    int   len;
} mfa_seq_t;

#define MFA_SEQ_BATCH 128

void mfa2sortedctxobj64_v2 (sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat)
{
    const bool resolve_conf = !sketch_opt_val->conflict;
    const int  nthreads     = sketch_opt_val->p;

    uint64_t totle_sketch_size = 0;   // cumulative across all sequences

    // ---- global index: one entry per "genome" (sequence), plus 0 at the front ----
    u64vec sketch_index;
    v_init(&sketch_index, 1024);      // initial capacity
    v_push(&sketch_index, 0ULL);      // index[0] = 0

    // tmpname: one name per genome (sequence)
    int tmpname_size_alloc = (int)sketch_index.cap;
    char (*tmpname)[PATHLEN] = malloc((size_t)tmpname_size_alloc * PATHLEN);
    if (!tmpname)
        err(errno, "%s(): Memory allocation failed for tmpname", __func__);

    // combined sketch file (keys only)
    FILE *comb_sketch_fp =
        fopen(format_string("%s/%s", sketch_opt_val->outdir, combined_sketch_suffix), "wb");
    if (!comb_sketch_fp)
        err(errno, "%s() open file error: %s/%s",
            __func__, sketch_opt_val->outdir, combined_sketch_suffix);
    setvbuf(comb_sketch_fp, NULL, _IOFBF, 8u << 20);

    // =========================
    //  loop over MFA *files*
    // =========================
    for (int fi = 0; fi < infile_stat->infile_num; ++fi) {
        char *seqfname = infile_stat->organized_infile_tab[fi].fpath;

        gzFile infile = gzopen(seqfname, "r");
        if (!infile)
            err(errno, "mfa2sortedctxobj64(): Cannot open file %s", seqfname);

        kseq_t *seq = kseq_init(infile);

        for (;;) {
            mfa_seq_t seqs[MFA_SEQ_BATCH];
            int       nseq = 0;
            int       ret;
            bool      eof = false;

            // -------- read up to MFA_SEQ_BATCH usable sequences into this batch --------
            while (nseq < MFA_SEQ_BATCH) {
                ret = kseq_read(seq);
                if (ret < 0) {  // EOF or error
                    eof = true;
                    break;
                }

                if (seq->seq.l <= klen)
                    continue;   // too short for any k-mer; don't treat as genome

                seqs[nseq].len = (int)seq->seq.l;

                // copy sequence
                seqs[nseq].seq = (char *)malloc((size_t)seq->seq.l);
                if (!seqs[nseq].seq)
                    err(errno, "%s(): OOM for seq string", __func__);
                memcpy(seqs[nseq].seq, seq->seq.s, (size_t)seq->seq.l);

                // copy name
                seqs[nseq].name = strdup(seq->name.s);
                if (!seqs[nseq].name)
                    err(errno, "%s(): OOM for seq name", __func__);

                ++nseq;
            }

            if (nseq == 0) {
                // no usable sequences in this batch
                if (eof) break;  // finished this file
                else continue;   // only short sequences encountered; keep reading
            }

            // -------- parallel sketch & KV build per sequence in this batch --------
            uint64_t **keys = (uint64_t **)calloc((size_t)nseq, sizeof(uint64_t *));
            uint64_t  *lens = (uint64_t  *)calloc((size_t)nseq, sizeof(uint64_t));
            if (!keys || !lens)
                err(errno, "%s(): OOM for per-sequence KV arrays", __func__);

            #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1)
            for (int si = 0; si < nseq; ++si) {
                u64vec vec;
                v_init(&vec, 1u << 15);

                sketch_read_into_vec(seqs[si].seq,seqs[si].len,&vec,
                                    ctxmask,tupmask,Bitslen.obj,klen);

                if (vec.n) {
                    // no abundance, no kmerocrs: just dedup + sort, get counts we discard
                    SortedKV_Arrays_t kv = build_kv_from_vec(&vec, false);

                    if (resolve_conf)
                        remove_ctx_with_conflict_obj_noabund(kv.keys, &kv.len, Bitslen.obj);

                    keys[si] = kv.keys;
                    lens[si] = kv.len;

                    free(kv.values);  // counts not kept in MFA mode
                } else {
                    lens[si] = 0;
                }

                v_free(&vec);
            }

            // -------- serial write in original sequence order for this batch --------
            for (int si = 0; si < nseq; ++si) {
                // ensure tmpname capacity (sketch_index.n == #genomes_so_far + 1)
                if ((int)sketch_index.n >= tmpname_size_alloc) {
                    tmpname_size_alloc += 1000;
                    char (*newtmp)[PATHLEN] =
                        realloc(tmpname, (size_t)tmpname_size_alloc * PATHLEN);
                    if (!newtmp)
                        err(errno, "%s(): Realloc failed for tmpname", __func__);
                    tmpname = newtmp;
                }

                // record this sequence name (one per genome)
                replace_special_chars_with_underscore(seqs[si].name);
                strncpy(tmpname[sketch_index.n - 1], seqs[si].name, PATHLEN);

                // write keys for this sequence
                if (lens[si])
                    fwrite(keys[si], sizeof(uint64_t), lens[si], comb_sketch_fp);

                totle_sketch_size += lens[si];
                v_push(&sketch_index, totle_sketch_size);

                free(keys[si]);
            }

            free(keys);
            free(lens);

            // free per-sequence buffers for this batch
            for (int si = 0; si < nseq; ++si) {
                free(seqs[si].seq);
                free(seqs[si].name);
            }

            if (eof) break;  // we've reached end of file
        } // end per-file batches

        kseq_destroy(seq);
        gzclose(infile);

        fprintf(stderr, "\r%dth/%d multifasta file %s completed!\t #genomes=%lu",
                fi + 1,
                infile_stat->infile_num,
                infile_stat->organized_infile_tab[fi].fpath,
                (unsigned long)(sketch_index.n - 1));
        if (fi == infile_stat->infile_num - 1) fprintf(stderr, "\n");
    } // for each file

    fclose(comb_sketch_fp);

    // write global index & stat
    write_to_file(test_create_fullpath(sketch_opt_val->outdir, idx_sketch_suffix),
                  sketch_index.a,sketch_index.n * sizeof(uint64_t));

    comblco_stat_one.infile_num = (uint32_t)(sketch_index.n - 1);
    concat_and_write_to_file(test_create_fullpath(sketch_opt_val->outdir, sketch_stat),
                             &comblco_stat_one,sizeof(comblco_stat_one),tmpname,
                             (size_t)comblco_stat_one.infile_num * PATHLEN);

    v_free(&sketch_index);
    free(tmpname);
}
