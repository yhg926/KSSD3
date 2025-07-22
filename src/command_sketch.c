#include "command_sketch.h"
#include "sketch_rearrange.h"
#include "kssdlib_sort.h"
#include <time.h>

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

uint32_t get_sketching_id(uint32_t hclen, uint32_t holen, uint32_t iolen, uint32_t drfold, uint32_t FILTER)
{
    return GET_SKETCHING_ID(hclen, holen, iolen, drfold, FILTER);
}

/*
void public_vars_init(dim_sketch_stat_t* sketch_stat_raw){
     FILTER = UINT32_MAX >> sketch_stat_raw->drfold  ;
     sketch_stat_raw->hash_id  =  GET_SKETCHING_ID(sketch_stat_raw->hclen, sketch_stat_raw->holen, sketch_stat_raw->iolen, sketch_stat_raw->drfold , FILTER);
     klen = 2*(sketch_stat_raw->hclen + sketch_stat_raw->holen) + sketch_stat_raw->iolen;
     if ( klen > 32 || FILTER < 256) err(EINVAL,"compute_sketch(): klen (%d) or FILTER (%u) is out of range (klen <=32 and FILTER: 256..0xffffffff)",klen, FILTER);
     comblco_stat_one = *sketch_stat_raw ;
     printf("Sketching method hashid = %u\tklen=%u\n", comblco_stat_one.hash_id,kcomblco_stat_one.len);
}
*/

void compute_sketch(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat)
{
    set_uint64kmer2generic_ctxobj(sketch_opt_val->coden_ctxobj_pattern);
    if (sketch_opt_val->split_mfa)
    { // mfa files parse
        mfa2sortedctxobj64(sketch_opt_val, infile_stat);
        return;
    }
    read_genomes2mem2sortedctxobj64(sketch_opt_val, infile_stat, 1000);
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

#include <x86intrin.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#ifndef MAP_POPULATE
#define MAP_POPULATE 0
#endif
#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS 0x20
#endif
#include <bits/mman-linux.h>

#define CACHE_ALIGN __attribute__((aligned(64)))
#define BATCH_SIZE (4 << 20) // 4MB批量处理
#define SIMD_WIDTH 32        // AVX2处理32字节块
#define KMER_BATCH 8         // 每个SIMD批次处理8个k-mer

// 缓存对齐的内存块结构
struct MemoryBlocks
{
    char seq_buffer[BATCH_SIZE] CACHE_ALIGN;
    uint64_t sketch_buffer[(BATCH_SIZE / 32) * KMER_BATCH] CACHE_ALIGN;
};

// SIMD辅助函数：将ASCII字符转换为2-bit编码
static inline __m256i simd_base_encode(__m256i input)
{
    // ASCII值处理：A(65)=00, C(67)=01, G(71)=10, T(84)=11
    const __m256i mask = _mm256_set1_epi8(0x1F);  // 取低5位
    const __m256i baseA = _mm256_set1_epi8(0x01); // A的掩码
    const __m256i baseC = _mm256_set1_epi8(0x03); // C的掩码
    const __m256i baseG = _mm256_set1_epi8(0x07); // G的掩码
    const __m256i baseT = _mm256_set1_epi8(0x14); // T的掩码

    __m256i masked = _mm256_and_si256(input, mask);
    __m256i cmpA = _mm256_cmpeq_epi8(masked, baseA);
    __m256i cmpC = _mm256_cmpeq_epi8(masked, baseC);
    __m256i cmpG = _mm256_cmpeq_epi8(masked, baseG);
    __m256i cmpT = _mm256_cmpeq_epi8(masked, baseT);

    __m256i result = _mm256_set1_epi8(0xFF); // 默认无效
    result = _mm256_blendv_epi8(result, _mm256_set1_epi8(0), cmpA);
    result = _mm256_blendv_epi8(result, _mm256_set1_epi8(1), cmpC);
    result = _mm256_blendv_epi8(result, _mm256_set1_epi8(2), cmpG);
    result = _mm256_blendv_epi8(result, _mm256_set1_epi8(3), cmpT);
    return result;
}

// 核心处理函数
int opt_seq2sortedsketch64(char *seqfname, char *outfname, bool abundance, int n)
{
    // 初始化内存池
    static struct MemoryBlocks *mem_pool = NULL;
    if (!mem_pool)
    {
        mem_pool = mmap(NULL, sizeof(struct MemoryBlocks),
                        PROT_READ | PROT_WRITE,
                        MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE, -1, 0);
    }

    // 优化gzip读取配置
    gzFile infile = gzopen(seqfname, "r");
    gzbuffer(infile, 1 << 20); // 1MB解压缓冲区
    kseq_t *seq = kseq_init(infile);

    // 预分配内存（调整为预估最大值的2倍）
    size_t max_sketch = (BATCH_SIZE / 32) * KMER_BATCH * 2; // 估算k-mer密度
    Vector raw_sketch;
    vector_init(&raw_sketch, sizeof(uint64_t));
    vector_reserve(&raw_sketch, max_sketch);

    uint64_t tuple = 0, crvstuple = 0;
    const uint32_t len_mv = 2 * klen - 2;
    const uint64_t ctxmask_local = ctxmask;
    const uint64_t tupmask_local = tupmask;

    // 批量读取处理循环
    while (kseq_read(seq) >= 0)
    {
        const char *s = seq->seq.s;
        const int seq_len = seq->seq.l;
        if (seq_len < klen)
            continue;

        int pos = 0;
        int base = 0;
#ifdef __AVX2__
        // SIMD处理32字节块
        for (; pos + SIMD_WIDTH <= seq_len; pos += SIMD_WIDTH)
        {
            __m256i chunk = _mm256_loadu_si256((__m256i *)(s + pos));
            __m256i bases = simd_base_encode(chunk);

            // 检测无效碱基
            __m256i invalid = _mm256_cmpeq_epi8(bases, _mm256_set1_epi8(0xFF));
            if (!_mm256_testz_si256(invalid, invalid))
            {
                base = 0; // 存在无效碱基则重置
                continue;
            }

            // 将向量转换为64位处理
            uint64_t *base64 = (uint64_t *)&bases;
            uint64_t current_tuple = tuple;
            uint64_t current_crvstuple = crvstuple;

            // 处理32字节中的每个碱基
            for (int i = 0; i < 4; ++i)
            { // 每个uint64_t包含8个碱基
                uint64_t block = base64[i];
                for (int j = 0; j < 8; ++j, block >>= 8)
                {
                    uint64_t b = block & 0xFF;
                    current_tuple = (current_tuple << 2) | b;
                    current_crvstuple = (current_crvstuple >> 2) | ((b ^ 3LLU) << len_mv);

                    if (++base >= klen)
                    {
                        // 生成k-mer特征
                        uint64_t mask_tuple = current_tuple & ctxmask_local;
                        uint64_t mask_crvstuple = current_crvstuple & ctxmask_local;
                        uint64_t unituple = (mask_tuple < mask_crvstuple) ? current_tuple : current_crvstuple;
                        unituple = uint64kmer2generic_ctxobj(unituple);

                        // 哈希过滤
                        uint64_t unictx = unituple & ctxmask_local;
                        if (SKETCH_HASH(unictx) <= FILTER)
                        {
                            if (raw_sketch.size >= max_sketch)
                            {
                                max_sketch *= 2;
                                vector_reserve(&raw_sketch, max_sketch);
                            }
                            uint64_t *dst = (uint64_t *)raw_sketch.data + raw_sketch.size++;
                            *dst = unituple & tupmask_local;
                        }
                    }
                }
            }
            tuple = current_tuple;
            crvstuple = current_crvstuple;
        }
#endif

        // 标量处理剩余部分
        for (; pos < seq_len; ++pos)
        {
            unsigned char c = s[pos];
            if (Basemap[c] == DEFAULT)
            {
                base = 0;
                tuple = 0;
                crvstuple = 0;
                continue;
            }
            uint64_t basenum = Basemap[c];
            uint64_t new_tuple = (tuple << 2) | basenum;
            uint64_t new_crvstuple = (crvstuple >> 2) | ((basenum ^ 3LLU) << len_mv);

            // 预取优化
            if ((pos & 0xF) == 0)
            {
                _mm_prefetch(s + pos + 64, _MM_HINT_T0);
            }

            if (++base >= klen)
            {
                uint64_t mask_tuple = new_tuple & ctxmask_local;
                uint64_t mask_crvstuple = new_crvstuple & ctxmask_local;
                uint64_t unictx, unituple;
                if (mask_tuple < mask_crvstuple)
                {
                    unictx = mask_tuple;
                    unituple = new_tuple;
                }
                else
                {
                    unictx = mask_crvstuple;
                    unituple = new_crvstuple;
                }
                // uint64_t unituple = (mask_tuple < mask_crvstuple) ? new_tuple : new_crvstuple;
                if (SKETCH_HASH(unictx) <= FILTER)
                {
                    if (raw_sketch.size >= max_sketch)
                    {
                        max_sketch *= 2;
                        vector_reserve(&raw_sketch, max_sketch);
                    }
                    uint64_t *dst = (uint64_t *)raw_sketch.data + raw_sketch.size++;
                    *dst = uint64kmer2generic_ctxobj(unituple);
                }
            }
            tuple = new_tuple;
            crvstuple = new_crvstuple;
        }
    }

    gzclose(infile);
    kseq_destroy(seq);

    // 并行排序优化
    uint64_t *mem_lco = raw_sketch.data;
    const size_t sketch_size = raw_sketch.size;
#pragma omp parallel
    {
#pragma omp single
        qsort(mem_lco, sketch_size, sizeof(uint64_t), qsort_comparator_uint64);
    }

    // 去重与过滤
    uint32_t *mem_ab = NULL;
    size_t kmer_ct = dedup_with_counts(mem_lco, sketch_size, &mem_ab);
    if (n > 1)
    {
        size_t write_pos = 0;
        for (size_t read_pos = 0; read_pos < kmer_ct; ++read_pos)
        {
            if (mem_ab[read_pos] >= n)
            {
                mem_lco[write_pos] = mem_lco[read_pos];
                mem_ab[write_pos] = mem_ab[read_pos];
                ++write_pos;
            }
        }
        kmer_ct = write_pos;
    }

    // 批量写入优化
    FILE *out = fopen(outfname, "wb");
    setvbuf(out, NULL, _IOFBF, 1 << 26); // 64MB写入缓冲区
    fwrite(mem_lco, sizeof(uint64_t), kmer_ct, out);
    fclose(out);

    if (abundance)
    {
        FILE *ab_file = fopen(format_string("%s.a", outfname), "wb");
        setvbuf(ab_file, NULL, _IOFBF, 1 << 26);
        fwrite(mem_ab, sizeof(uint32_t), kmer_ct, ab_file);
        fclose(ab_file);
        free(mem_ab);
    }
    vector_free(&raw_sketch);
    return kmer_ct;
}

// 新增预定义参数
#define CACHE_LINE_SIZE 64
#define GZ_BUFFER_SIZE (1 << 20)   // 1MB解压缓冲区
#define INIT_SKETCH_SIZE (1 << 24) // 初始预分配16M元素
#define BATCH_PROCESS_SIZE 1024    // 批量处理单位

int opt2_seq2sortedsketch64(char *seqfname, char *outfname, bool abundance, int n)
{ // no SIMD optimization
    // 优化1：设置更大的gzip解压缓冲区
    gzFile infile = gzopen(seqfname, "r");
    if (!infile)
        err(errno, "reads2sketch64(): Cannot open file %s", seqfname);
    gzbuffer(infile, GZ_BUFFER_SIZE);

    kseq_t *seq = kseq_init(infile);

    // 优化2：预分配对齐内存
    uint64_t *raw_sketch = aligned_alloc(CACHE_LINE_SIZE, INIT_SKETCH_SIZE * sizeof(uint64_t));
    size_t sketch_capacity = INIT_SKETCH_SIZE;
    size_t sketch_size = 0;

    // 优化3：批量处理中间变量
    uint64_t batch_sketch[BATCH_PROCESS_SIZE] __attribute__((aligned(CACHE_LINE_SIZE)));
    size_t batch_count = 0;

    uint64_t tuple = 0, crvstuple = 0;
    const uint32_t len_mv = 2 * klen - 2;

    // 优化4：提前计算常用掩码
    const uint64_t ctxmask_local = ctxmask;
    const uint64_t tupmask_local = tupmask;

    while (kseq_read(seq) >= 0)
    {
        const char *s = seq->seq.s;
        const int seq_len = seq->seq.l;
        if (seq_len < klen)
            continue;

        // 优化5：循环展开与局部变量缓存
        int base = 0;
        uint64_t current_tuple = tuple;
        uint64_t current_crvstuple = crvstuple;

        for (int pos = 0; pos < seq_len; pos++)
        {
            const unsigned char c = (unsigned short)s[pos];
            const uint64_t basenum = Basemap[c];

            if (basenum == DEFAULT)
            {
                base = 0;
                current_tuple = 0;
                current_crvstuple = 0;
                continue;
            }

            // 优化6：拆解依赖链
            const uint64_t new_tuple = (current_tuple << 2) | basenum;
            const uint64_t new_crvstuple = (current_crvstuple >> 2) |
                                           ((basenum ^ 3LLU) << len_mv);

            if (++base >= klen)
            {
                // 优化7：无分支条件选择
                const uint64_t masked_tuple = new_tuple & ctxmask_local;
                const uint64_t masked_crvstuple = new_crvstuple & ctxmask_local;
                const uint64_t unituple = (masked_tuple < masked_crvstuple) ? new_tuple : new_crvstuple;

                // 优化8：批量写入
                batch_sketch[batch_count++] = uint64kmer2generic_ctxobj(unituple) & tupmask_local;

                // 批量提交
                if (batch_count == BATCH_PROCESS_SIZE)
                {
                    if (sketch_size + batch_count > sketch_capacity)
                    {
                        // 优化9：指数扩容策略
                        sketch_capacity = (sketch_size + batch_count) * 2;
                        raw_sketch = realloc(raw_sketch, sketch_capacity * sizeof(uint64_t));
                    }
                    memcpy(raw_sketch + sketch_size, batch_sketch,
                           batch_count * sizeof(uint64_t));
                    sketch_size += batch_count;
                    batch_count = 0;
                }
            }

            current_tuple = new_tuple;
            current_crvstuple = new_crvstuple;
        }

        // 更新全局变量供下次循环使用
        tuple = current_tuple;
        crvstuple = current_crvstuple;
    }
    printf("OK0\n");
    // 处理剩余批次
    if (batch_count > 0)
    {
        if (sketch_size + batch_count > sketch_capacity)
        {
            sketch_capacity = sketch_size + batch_count;
            raw_sketch = realloc(raw_sketch, sketch_capacity * sizeof(uint64_t));
        }
        memcpy(raw_sketch + sketch_size, batch_sketch,
               batch_count * sizeof(uint64_t));
        sketch_size += batch_count;
    }
    printf("OK1\n");
    gzclose(infile);
    kseq_destroy(seq);

    // 排序与去重
    qsort(raw_sketch, sketch_size, sizeof(uint64_t), qsort_comparator_uint64);
    uint32_t *mem_ab = NULL;
    const uint32_t unique_count = dedup_with_counts(raw_sketch, sketch_size, &mem_ab);
    printf("OK2\n");
    // 过滤低频k-mer
    uint32_t final_count = unique_count;
    if (n > 1)
    {
        final_count = 0;
        for (uint32_t i = 0; i < unique_count; i++)
        {
            if (mem_ab[i] >= n)
            {
                raw_sketch[final_count] = raw_sketch[i];
                mem_ab[final_count] = mem_ab[i];
                final_count++;
            }
        }
    }
    printf("OK3\n");
    // 输出处理
    write_to_file(outfname, raw_sketch, final_count * sizeof(uint64_t));
    if (abundance)
    {
        write_to_file(format_string("%s.a", outfname), mem_ab, final_count * sizeof(uint32_t));
        free(mem_ab);
    }

    free(raw_sketch);
    return final_count;
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
                if (len < klen) continue;
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
                    if (++base < klen) continue;

                    unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
                    unictx = unituple & ctxmask;
                    //printf("%lx\t%lx\t%lx\t%lx\t%lx\t%lx\n", unituple, unictx, tuple,tuple & ctxmask, crvstuple, crvstuple & ctxmask);
                    if (SKETCH_HASH(unictx) > FILTER)
                        continue;

                    int ret;
                    khint_t key = kh_put(sort64, h, uint64kmer2generic_ctxobj(unituple & tupmask), &ret);
                    if (ret)
                        kh_value(h, key) = 1;
                } // nt pos loop
                free(s);
            } // seq j loop
            SortedKV_Arrays_t lco_ab = sort_khash_u64(h);
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

        printf("\r%d/%d genome proceed", infile_num_p, infile_stat->infile_num);

    } // infile_num_p loop

    for (int i = 0; i < infile_stat->infile_num; i++)
        sketch_index[i + 1] += sketch_index[i];
    write_to_file(format_string("%s/%s", sketch_opt_val->outdir, idx_sketch_suffix), sketch_index, (infile_stat->infile_num + 1) * sizeof(sketch_index[0]));

    fclose(comb_sketch_fp);
    vector_free(&all_reads);
    free_all(sketch_index, gseq_nums, batch_sketches, NULL);
    write_sketch_stat(sketch_opt_val->outdir, infile_stat);
}
