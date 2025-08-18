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
void test_read_genomes2mem2sortedctxobj64(sketch_opt_t *sketch_opt_val,
                                          infile_tab_t *infile_stat,
                                          int batch_size);
uint32_t get_sketching_id(uint32_t hclen, uint32_t holen, uint32_t iolen, uint32_t drfold, uint32_t FILTER)
{
    return GET_SKETCHING_ID(hclen, holen, iolen, drfold, FILTER);
}

void compute_sketch(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat)
{
    if (sketch_opt_val->split_mfa)
    { // mfa files parse
        mfa2sortedctxobj64(sketch_opt_val, infile_stat);
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
simple_sketch_t *simple_genomes2mem2sortedctxobj64_mem(infile_tab_t *infile_stat, int drfold)
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

static inline void apply_conflict_filter(SortedKV_Arrays_t *kv,bool has_abundance,
                                         bool conflict_flag,uint32_t nobjbits)
{
    if (conflict_flag || kv->len == 0) return;
    if (has_abundance) {
        remove_ctx_with_conflict_obj(kv, nobjbits);
    } else {
        remove_ctx_with_conflict_obj_noabund(kv->keys, &kv->len, nobjbits);
    }
}

static inline void write_index_payload_and_free(FILE *comb,FILE *comb_ab,
                bool has_abundance,uint64_t *sketch_index_slot, SortedKV_Arrays_t *kv)
{
    *sketch_index_slot = kv->len;
    if (kv->len) {
        fwrite(kv->keys,   sizeof(uint64_t), kv->len, comb);
        if (has_abundance) fwrite(kv->values, sizeof(uint32_t), kv->len, comb_ab);
    }
    free(kv->keys);   kv->keys   = NULL;
    free(kv->values); kv->values = NULL;
    kv->len = 0;
}

// ---------------- Mode A: per-genome parallel (vector + radix + counts) ----------------
static void sketch_many_files_in_parallel(sketch_opt_t *opt, infile_tab_t *tab, int batch_size)
{
    const int nfiles = tab->infile_num;
    const bool has_abundance = opt->abundance;

    FILE *comb = fopen(format_string("%s/%s", opt->outdir, combined_sketch_suffix), "wb");
    if (!comb) err(errno, "%s() open comb", __func__);
    setvbuf(comb, NULL, _IOFBF, 8u << 20);

    FILE *comb_ab = NULL;
    if (has_abundance) {
        comb_ab = fopen(format_string("%s/%s", opt->outdir, combined_ab_suffix), "wb");
        if (!comb_ab) err(errno, "%s() open file error: %s/%s", __func__, opt->outdir, combined_ab_suffix);
        setvbuf(comb_ab, NULL, _IOFBF, 8u << 20);
    }

    uint64_t *sketch_index = (uint64_t *)calloc((size_t)nfiles + 1, sizeof(uint64_t));
    if (!sketch_index) err(errno, "%s(): OOM index", __func__);

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

        #pragma omp parallel for num_threads(opt->p) schedule(dynamic, 1)
        for (int bi = 0; bi < this_batch; ++bi) {
            const int file_idx = batch_start + bi;
            const char *fpath  = tab->organized_infile_tab[file_idx].fpath;

            gzFile infile = gzopen(fpath, "r");
            if (!infile) err(errno, "%s(): Cannot open %s", __func__, fpath);
            (void)gzbuffer(infile, 4u << 20);

            kseq_t *seq = kseq_init(infile);
            if (!seq) err(errno, "%s(): kseq_init %s", __func__, fpath);

            u64vec vec; v_init(&vec, 1u << 15);
            while (kseq_read(seq) >= 0) {
                sketch_read_into_vec(seq->seq.s, (int)seq->seq.l,
                                     &vec, ctxmask, tupmask, Bitslen.obj, klen);
            }

            if (vec.n) {
                SortedKV_Arrays_t kv = build_kv_from_vec(&vec, has_abundance);
                apply_conflict_filter(&kv, has_abundance, opt->conflict, Bitslen.obj);
                batch_keys[bi] = kv.keys;
                lens[bi]       = kv.len;
                if (has_abundance) batch_vals[bi] = kv.values;
                else               free(kv.values);
            } else {
                lens[bi] = 0;
            }

            v_free(&vec);
            kseq_destroy(seq);
            gzclose(infile);
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

        free_all(batch_keys,lens,NULL);
        if (has_abundance) free(batch_vals);

        fprintf(stderr, "\r%d/%d genomes processed", batch_end, nfiles);
        fflush(stderr);
    }

    for (int i = 0; i < nfiles; ++i) sketch_index[i + 1] += sketch_index[i];
    write_to_file(format_string("%s/%s", opt->outdir, idx_sketch_suffix),
                  sketch_index, (size_t)(nfiles + 1) * sizeof(uint64_t));
    fclose(comb);
    if (has_abundance) fclose(comb_ab);
    free(sketch_index);
    write_sketch_stat(opt->outdir, tab);
}

// ---------------- Mode B: vector + bucketed parallel finalizer (with counts) ----------------
static inline unsigned clamp_bucket_bits_from_total(size_t total_distinct){
    // target ~4k entries per bucket; clamp BITS to [8..14]
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

static void sketch_few_files_with_intrafile_parallel(sketch_opt_t *opt, infile_tab_t *tab, int BATCH_READS)
{
    const int nfiles = tab->infile_num;
    const bool has_abundance = true; // Mode-B keeps counts; we can still write only keys if desired

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
        gzFile in = gzopen(path, "r");
        if (!in) err(errno, "%s(): Cannot open %s", __func__, path);
        (void)gzbuffer(in, 4u<<20);
        kseq_t *seq = kseq_init(in);

        typedef struct { char *s; int l; } read_span_t;
        read_span_t *R = (read_span_t*)malloc((size_t)BATCH_READS*sizeof(*R));
        if (!R) err(errno, "%s(): OOM batch reads", __func__);
        int rcnt=0;

        const int nth = (opt->p > 0 ? opt->p : 1);
        u64vec *V_thr = (u64vec*)malloc((size_t)nth*sizeof(u64vec));
        if (!V_thr) err(errno, "%s(): OOM V_thr", __func__);
        for (int t=0;t<nth;++t) v_init(&V_thr[t], 1u<<15);

        uint32_t **C_thr = (uint32_t**)calloc((size_t)nth, sizeof(uint32_t*));
        if (!C_thr) err(errno, "%s(): OOM C_thr", __func__);

        // read → per-thread vectors
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
                #pragma omp parallel for if (nth>1) num_threads(opt->p) schedule(static)
                for (int ri=0; ri<rcnt; ++ri){
                    int tid = 0;
                    #ifdef _OPENMP
                    tid = omp_get_thread_num();
                    #endif
                    sketch_read_into_vec(R[ri].s, R[ri].l, &V_thr[tid],
                                         ctxmask, tupmask, Bitslen.obj, klen);
                }
                for (int i=0;i<rcnt;++i) free(R[i].s);
                rcnt=0;
            }
        } while (ret >= 0);
        free(R);

        // per-thread shrink: sort + counts
        size_t pre_elems=0; for (int t=0;t<nth;++t) pre_elems += V_thr[t].n;
        #pragma omp parallel for if (nth>1 && pre_elems >= (1u<<20)) num_threads(opt->p) schedule(dynamic)
        for (int t=0; t<nth; ++t) shrink_thread_vec(&V_thr[t], &C_thr[t]);

        size_t total_distinct=0; for (int t=0;t<nth;++t) total_distinct += V_thr[t].n;
        if (total_distinct==0){
            sketch_index[f+1]=0;
            for (int t=0;t<nth;++t){ v_free(&V_thr[t]); free(C_thr[t]); }
            free(V_thr); free(C_thr);
            kseq_destroy(seq); gzclose(in);
            fprintf(stderr, "\r%d/%d genomes processed", f+1, nfiles); fflush(stderr);
            continue;
        }

        // buckets
        const unsigned BITS = clamp_bucket_bits_from_total(total_distinct);
        const size_t NB = (size_t)1u << BITS;

        size_t *bucket_totals = (size_t*)calloc(NB, sizeof(size_t));
        size_t **thr_counts   = (size_t**)malloc((size_t)nth * sizeof(size_t*));
        if (!bucket_totals || !thr_counts) err(errno, "%s(): OOM bucket meta", __func__);
        for (int t=0;t<nth;++t){
            thr_counts[t] = (size_t*)calloc(NB, sizeof(size_t));
            if (!thr_counts[t]) err(errno, "%s(): OOM thr_counts", __func__);
        }

        #pragma omp parallel for num_threads(opt->p) schedule(static)
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
        if (!bucket_off) err(errno, "%s(): OOM bucket_off", __func__);
        bucket_off[0]=0;
        for (size_t b=0;b<NB;++b) bucket_off[b+1] = bucket_off[b] + bucket_totals[b];

        const size_t TOTAL_ENTRIES = bucket_off[NB];
        uint64_t *B_keys = (uint64_t*)malloc(TOTAL_ENTRIES*sizeof(uint64_t));
        uint32_t *B_vals = (uint32_t*)malloc(TOTAL_ENTRIES*sizeof(uint32_t));
        if (!B_keys || !B_vals) err(errno, "%s(): OOM bucket KV", __func__);

        size_t **thr_write = (size_t**)malloc((size_t)nth * sizeof(size_t*));
        if (!thr_write) err(errno, "%s(): OOM thr_write", __func__);
        for (int t=0;t<nth;++t){
            thr_write[t] = (size_t*)malloc(NB*sizeof(size_t));
            if (!thr_write[t]) err(errno, "%s(): OOM thr_write[t]", __func__);
        }
        for (size_t b=0;b<NB;++b){
            size_t off = bucket_off[b];
            for (int t=0;t<nth;++t){ thr_write[t][b] = off; off += thr_counts[t][b]; }
        }

        #pragma omp parallel for num_threads(opt->p) schedule(static)
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

        // per-bucket parallel KV sort+merge
        uint64_t **bk_keys = (uint64_t**)malloc(NB*sizeof(uint64_t*));
        uint32_t **bk_vals = (uint32_t**)malloc(NB*sizeof(uint32_t*));
        size_t    *bk_len  = (size_t*)   malloc(NB*sizeof(size_t));
        if (!bk_keys || !bk_vals || !bk_len) err(errno, "%s(): OOM bucket outs", __func__);

        #pragma omp parallel for num_threads(opt->p) schedule(dynamic)
        for (size_t b=0; b<NB; ++b){
            const size_t begin = bucket_off[b], end = bucket_off[b+1];
            const size_t n = end - begin;
            if (!n){ bk_keys[b]=NULL; bk_vals[b]=NULL; bk_len[b]=0; continue; }
            uint64_t *kseg = B_keys + begin; uint32_t *vseg = B_vals + begin;
            radix_sort_kv_u64(kseg, vseg, n);
            size_t m = shrink_kv_inplace_u64_u32(kseg, vseg, n);
            uint64_t *ko = (uint64_t*)malloc(m*sizeof(uint64_t));
            uint32_t *vo = (uint32_t*)malloc(m*sizeof(uint32_t));
            if (!ko || !vo) err(errno, "%s(): OOM bucket kv out", __func__);
            memcpy(ko, kseg, m*sizeof(uint64_t));
            memcpy(vo, vseg, m*sizeof(uint32_t));
            bk_keys[b]=ko; bk_vals[b]=vo; bk_len[b]=m;
        }

        free(B_keys); free(B_vals); free(bucket_totals);

        size_t *out_off = (size_t*)malloc((NB+1)*sizeof(size_t));
        if (!out_off) err(errno, "%s(): OOM out_off", __func__);
        out_off[0]=0; for (size_t b=0;b<NB;++b) out_off[b+1] = out_off[b] + bk_len[b];
        const size_t M = out_off[NB];

        uint64_t *out_keys = (uint64_t*)malloc(M*sizeof(uint64_t));
        uint32_t *out_vals = (uint32_t*)malloc(M*sizeof(uint32_t));
        if (!out_keys || !out_vals) err(errno, "%s(): OOM final kv", __func__);

        #pragma omp parallel for num_threads(opt->p) schedule(static)
        for (size_t b=0; b<NB; ++b){
            if (!bk_len[b]) continue;
            memcpy(out_keys + out_off[b], bk_keys[b], bk_len[b]*sizeof(uint64_t));
            memcpy(out_vals + out_off[b], bk_vals[b], bk_len[b]*sizeof(uint32_t));
            free(bk_keys[b]); free(bk_vals[b]);
        }
        free(bk_keys); free(bk_vals); free(bk_len); free(bucket_off); free(out_off);

        SortedKV_Arrays_t kv = (SortedKV_Arrays_t){ .keys=out_keys, .values=out_vals, .len=M };
        apply_conflict_filter(&kv, /*has_abundance=*/true, opt->conflict, Bitslen.obj);
        write_index_payload_and_free(comb, (opt->abundance?comb_ab:NULL), opt->abundance,
                                     &sketch_index[f+1], &kv);

        kseq_destroy(seq); gzclose(in);
        fprintf(stderr, "\r%d/%d genomes processed", f+1, nfiles); fflush(stderr);
    }

    for (int i=0;i<nfiles;++i) sketch_index[i+1] += sketch_index[i];
    write_to_file(format_string("%s/%s", opt->outdir, idx_sketch_suffix),
                  sketch_index, (size_t)(nfiles + 1) * sizeof(uint64_t));
    fclose(comb);
    if (opt->abundance) fclose(comb_ab);
    free(sketch_index);
    write_sketch_stat(opt->outdir, tab);
}

