// kssd3_sketch_fast.c
// High-throughput, rolling-hash FracMinHash sketcher core.
// Build: gcc -O3 -march=native -flto -pthread -DZLIB -lz -o ... (or -DLIBDEFLATE + -ldeflate)

#define _GNU_SOURCE
#include "kssd3_sketch_fast.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>
#include <inttypes.h>
#include <unistd.h>

// --- Optional gzip backends --------------------------------------------------
#if defined(LIBDEFLATE)
  #include <libdeflate.h>
#elif defined(ZLIB)
  #include <zlib.h>
#endif

// --- Utility macros ----------------------------------------------------------
#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

#define DIEF(fmt, ...) do { fprintf(stderr, "[sketch_fast] " fmt "\n", __VA_ARGS__); exit(2);} while(0)
#define DIE(msg)     do { fprintf(stderr, "[sketch_fast] %s\n", msg); exit(2);} while(0)

static inline int is_pow2_u32(uint32_t x) { return x && !(x & (x - 1)); }

// --- ASCII -> 2-bit LUT ------------------------------------------------------
// A,C,G,T -> 0,1,2,3 ; anything else = 0xFF (invalid)
static uint8_t DNA_LUT[256];
__attribute__((constructor))
static void init_lut(void) {
    memset(DNA_LUT, 0xFF, sizeof(DNA_LUT));
    DNA_LUT[(unsigned)'A'] = DNA_LUT[(unsigned)'a'] = 0;
    DNA_LUT[(unsigned)'C'] = DNA_LUT[(unsigned)'c'] = 1;
    DNA_LUT[(unsigned)'G'] = DNA_LUT[(unsigned)'g'] = 2;
    DNA_LUT[(unsigned)'T'] = DNA_LUT[(unsigned)'t'] = 3;
}

// --- Rolling 2-bit k-mer state (forward & reverse-complement) ----------------
// Supports k up to 31 (2*k bits <= 62); adjust to 63-bit if needed.
typedef struct {
    uint64_t fwd;      // 2-bit packed forward k-mer (lowest 2 bits = most recent base)
    uint64_t rev;      // 2-bit packed reverse-complement
    int      len;      // current length (<= k)
} roll2_t;

static inline uint64_t rc2(uint8_t b) { return (uint64_t)(b ^ 0x3u); } // A<->T, C<->G

static inline void roll2_reset(roll2_t *r) { r->fwd = r->rev = 0; r->len = 0; }

static inline void roll2_push(roll2_t *r, uint8_t b, int k) {
    // push base b (0..3) into forward; RC updates from the other side
    r->fwd = ((r->fwd << 2) | b) & ((k>=32)?~0ULL:((1ULL << (2*k)) - 1ULL));
    // For reverse-complement: shift right and insert complemented base at top
    r->rev = (r->rev >> 2) | (rc2(b) << (2*(k-1)));
    if (r->len < k) r->len++;
}

static inline uint64_t roll2_canonical(const roll2_t *r) {
    return (r->fwd < r->rev) ? r->fwd : r->rev;
}

// --- Hashing: mix 2-bit packed canonical k-mer to 64-bit ---------------------
// SplitMix64-style mixer (fast, well distributed)
static inline uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// --- FracMinHash keep test ---------------------------------------------------
typedef struct {
    uint32_t c;
    uint32_t mask;     // c-1 if power of two, else 0
    uint64_t thresh;   // threshold for non-power-of-two: keep if h < thresh
    int      use_mask;
} fmh_t;

static inline void fmh_init(fmh_t *f, uint32_t c) {
    f->c = c;
    if (is_pow2_u32(c)) {
        f->use_mask = 1; f->mask = c - 1; f->thresh = 0;
    } else {
        f->use_mask = 0; f->mask = 0;
        // We want P(keep)=1/c: keep if h < floor(2^64 / c)
        __uint128_t one = (((__uint128_t)1) << 64);
        f->thresh = (uint64_t)(one / c);
    }
}

static inline int fmh_keep(const fmh_t *f, uint64_t h) {
    return f->use_mask ? ((h & f->mask) == 0) : (h < f->thresh);
}

// --- Flat hash table (open addressing, linear probe, robin-hood-ish) ---------
typedef struct {
    uint64_t *keys;
    uint32_t *vals;
    uint8_t  *used;     // 0 empty, 1 used
    size_t    cap;
    size_t    n;
} map_t;

static map_t map_make(size_t cap) {
    map_t m = {0};
    m.cap = cap ? cap : 1024;
    m.keys = (uint64_t*)calloc(m.cap, sizeof(uint64_t));
    m.vals = (uint32_t*)calloc(m.cap, sizeof(uint32_t));
    m.used = (uint8_t*) calloc(m.cap, sizeof(uint8_t));
    if (!m.keys || !m.vals || !m.used) DIE("OOM map_make");
    return m;
}

static void map_free(map_t *m) {
    free(m->keys); free(m->vals); free(m->used);
    m->keys = NULL; m->vals = NULL; m->used = NULL; m->cap = m->n = 0;
}

static inline size_t map_probe(const map_t *m, uint64_t k) {
    // 64-bit Fibonacci hash to index
    size_t i = (size_t)((k * 11400714819323198485ull) >> (64 - 16)); // 2^16 buckets step; resized anyway
    if (unlikely(i >= m->cap)) i %= m->cap;
    return i;
}

static inline void map_grow(map_t *m) {
    size_t newcap = m->cap * 2;
    uint64_t *ok = m->keys; uint32_t *ov = m->vals; uint8_t *ou = m->used; size_t oc = m->cap;
    *m = map_make(newcap);
    for (size_t i=0;i<oc;i++) if (ou[i]) {
        uint64_t k = ok[i]; uint32_t v = ov[i];
        size_t p = map_probe(m, k);
        while (m->used[p]) { p = (p + 1) & (m->cap - 1); } // cap is power-of-two implicitly after doubling
        m->used[p] = 1; m->keys[p] = k; m->vals[p] = v; m->n++;
    }
    free(ok); free(ov); free(ou);
}

static inline void map_inc(map_t *m, uint64_t k) {
    if (unlikely((m->n * 2) >= m->cap)) map_grow(m);
    size_t p = map_probe(m, k);
    while (1) {
        if (!m->used[p]) { m->used[p]=1; m->keys[p]=k; m->vals[p]=1; m->n++; return; }
        if (m->keys[p]==k) { m->vals[p]++; return; }
        p = (p + 1) & (m->cap - 1);
    }
}

// --- Sharding ---------------------------------------------------------------
typedef struct {
    map_t *maps;        // array of maps, length = shards_per_thread
    size_t shards;
} shardset_t;

// --- Ring buffer for producer/consumer --------------------------------------
typedef struct {
    uint8_t *buf;
    size_t   cap;
    size_t   sz;
} chunk_t;

typedef struct {
    chunk_t  *chunks;
    size_t    n_chunks;
    size_t    head, tail;
    int       closed;
    pthread_mutex_t mu;
    pthread_cond_t  cv_put, cv_get;
} ring_t;

static void ring_init(ring_t *r, size_t n_chunks, size_t chunk_mb) {
    r->n_chunks = n_chunks; r->head = r->tail = 0; r->closed = 0;
    r->chunks = (chunk_t*)calloc(n_chunks, sizeof(chunk_t));
    if (!r->chunks) DIE("OOM ring_init");
    for (size_t i=0;i<n_chunks;i++) {
        r->chunks[i].cap = chunk_mb * 1024u * 1024u;
        r->chunks[i].buf = (uint8_t*)malloc(r->chunks[i].cap);
        if (!r->chunks[i].buf) DIE("OOM ring chunk");
        r->chunks[i].sz = 0;
    }
    pthread_mutex_init(&r->mu, NULL);
    pthread_cond_init(&r->cv_put, NULL);
    pthread_cond_init(&r->cv_get, NULL);
}

static void ring_free(ring_t *r) {
    for (size_t i=0;i<r->n_chunks;i++) free(r->chunks[i].buf);
    free(r->chunks);
    pthread_mutex_destroy(&r->mu);
    pthread_cond_destroy(&r->cv_put);
    pthread_cond_destroy(&r->cv_get);
}

static chunk_t* ring_get_putslot(ring_t *r) {
    pthread_mutex_lock(&r->mu);
    while (((r->head + 1) % r->n_chunks) == r->tail) pthread_cond_wait(&r->cv_put, &r->mu);
    chunk_t *c = &r->chunks[r->head];
    c->sz = 0;
    pthread_mutex_unlock(&r->mu);
    return c;
}

static void ring_commit_put(ring_t *r) {
    pthread_mutex_lock(&r->mu);
    r->head = (r->head + 1) % r->n_chunks;
    pthread_cond_signal(&r->cv_get);
    pthread_mutex_unlock(&r->mu);
}

static chunk_t* ring_get_getslot(ring_t *r) {
    pthread_mutex_lock(&r->mu);
    while ((r->tail == r->head) && !r->closed) pthread_cond_wait(&r->cv_get, &r->mu);
    if ((r->tail == r->head) && r->closed) { pthread_mutex_unlock(&r->mu); return NULL; }
    chunk_t *c = &r->chunks[r->tail];
    pthread_mutex_unlock(&r->mu);
    return c;
}

static void ring_release_get(ring_t *r) {
    pthread_mutex_lock(&r->mu);
    r->tail = (r->tail + 1) % r->n_chunks;
    pthread_cond_signal(&r->cv_put);
    pthread_mutex_unlock(&r->mu);
}

static void ring_close(ring_t *r) {
    pthread_mutex_lock(&r->mu);
    r->closed = 1;
    pthread_cond_broadcast(&r->cv_get);
    pthread_mutex_unlock(&r->mu);
}

// --- Minimal FASTA/FASTQ scanner over bytes ---------------------------------
// We decode bases and slide a k window across each record; invalid bases reset window.
typedef struct worker_ctx_t {
    const kssd3_sketch_cfg *cfg;
    fmh_t       fmh;
    shardset_t  ss;            // per-thread shard maps
    size_t      shard_bits;
    ring_t     *ring;
} worker_ctx_t;

static inline int guess_is_fastq(const uint8_t *buf, size_t n) {
    // crude: fastq starts with '@' and contains '+' on a later line; otherwise assume fasta '>'
    for (size_t i=0;i<n && i<4096;i++) {
        if (buf[i]=='>') return 0;
        if (buf[i]=='@') return 1;
    }
    return 1;
}

static void* worker_main(void *arg) {
    worker_ctx_t *w = (worker_ctx_t*)arg;
    const int k = w->cfg->k;
    const size_t shard_bits = w->shard_bits;
    const size_t nshards = (size_t)1 << shard_bits;

    // init shard maps
    w->ss.shards = nshards;
    w->ss.maps = (map_t*)malloc(nshards * sizeof(map_t));
    if (!w->ss.maps) DIE("OOM shard maps");
    // Rough capacity estimate: assume 1/c kept, reserve low LF per shard
    size_t per_shard_cap = 4096; // auto-grow if needed
    for (size_t i=0;i<nshards;i++) w->ss.maps[i] = map_make(per_shard_cap);

    while (1) {
        chunk_t *c = ring_get_getslot(w->ring);
        if (!c) break;
        const uint8_t *p = c->buf, *end = c->buf + c->sz;

        int mode_fastq = w->cfg->is_fastq >= 0 ? w->cfg->is_fastq : guess_is_fastq(p, (size_t)(end - p));

        roll2_t r; roll2_reset(&r);
        int in_seq = 0;
        size_t line_start = 0, i = 0, line_no = 0, fastq_phase = 0;

        while (i < (size_t)(end - c->buf)) {
            uint8_t ch = p[i++];
            if (mode_fastq) {
                // FASTQ: phases 0=@header,1=seq,2=+,3=qual
                if (ch == '\n') {
                    line_no++;
                    if (fastq_phase == 0) { fastq_phase = 1; in_seq = 1; roll2_reset(&r); }
                    else if (fastq_phase == 1) { fastq_phase = 2; in_seq = 0; }
                    else if (fastq_phase == 2) { fastq_phase = 3; in_seq = 0; }
                    else if (fastq_phase == 3) { fastq_phase = 0; in_seq = 0; }
                    line_start = i;
                    continue;
                }
                if (fastq_phase == 1 && ch != '\r') {
                    uint8_t b = DNA_LUT[ch];
                    if (b > 3) { roll2_reset(&r); continue; }
                    roll2_push(&r, b, k);
                    if (r.len == k) {
                        uint64_t key = mix64(roll2_canonical(&r));
                        if (fmh_keep(&w->fmh, key)) {
                            // shard by top bits
                            size_t shard = (size_t)(key >> (64 - shard_bits));
                            map_inc(&w->ss.maps[shard], key);
                        }
                    }
                }
            } else {
                // FASTA: '>' header lines; sequences in other lines
                if (ch == '>') { in_seq = 0; roll2_reset(&r); }
                else if (ch == '\n') { in_seq = 1; }
                else if (in_seq && ch != '\r') {
                    uint8_t b = DNA_LUT[ch];
                    if (b > 3) { roll2_reset(&r); continue; }
                    roll2_push(&r, b, k);
                    if (r.len == k) {
                        uint64_t key = mix64(roll2_canonical(&r));
                        if (fmh_keep(&w->fmh, key)) {
                            size_t shard = (size_t)(key >> (64 - shard_bits));
                            map_inc(&w->ss.maps[shard], key);
                        }
                    }
                }
            }
        }

        ring_release_get(w->ring);
    }

    return NULL;
}

// Emit all shards of all threads
typedef struct {
    size_t      nthreads;
    worker_ctx_t **workers;
    kssd3_emit_fn emit;
    void       *emit_user;
} emit_ctx_t;

static void emit_all(const emit_ctx_t *e) {
    for (size_t t=0; t<e->nthreads; ++t) {
        shardset_t *ss = &e->workers[t]->ss;
        for (size_t s=0; s<ss->shards; ++s) {
            map_t *m = &ss->maps[s];
            for (size_t i=0;i<m->cap;i++) if (m->used[i]) {
                e->emit(m->keys[i], m->vals[i], e->emit_user);
            }
        }
    }
}

// --- I/O thread: reads (plain or gz) into ring --------------------------------
typedef struct {
    ring_t *ring;
    const kssd3_sketch_cfg *cfg;
    int err;
} io_ctx_t;

static void* io_main(void *arg) {
    io_ctx_t *io = (io_ctx_t*)arg;
    const char *path = io->cfg->input_path;
    size_t chunk_mb = io->cfg->io_buf_mb ? io->cfg->io_buf_mb : 4;

    FILE *fp = NULL;
    int is_stdin = (strcmp(path, "-")==0);

#if defined(ZLIB)
    gzFile gzf = NULL;
#elif defined(LIBDEFLATE)
    // libdeflate raw stream: for simplicity, we fallback to zlib or plain in this sample
#endif

    if (is_stdin) {
#if defined(ZLIB)
        gzf = gzdopen(fileno(stdin), "rb");
        if (!gzf) { io->err=1; return NULL; }
#else
        fp = stdin;
#endif
    } else {
        // sniff .gz extension
        size_t len = strlen(path);
        int is_gz = (len>=3 && path[len-3]=='.' && path[len-2]=='g' && path[len-1]=='z');
#if defined(ZLIB)
        if (is_gz) {
            gzf = gzopen(path, "rb");
            if (!gzf) { io->err=1; return NULL; }
            gzbuffer(gzf, (unsigned)(chunk_mb*1024*1024));
        } else
#endif
        {
            fp = fopen(path, "rb");
            if (!fp) { io->err=1; return NULL; }
        }
    }

    for (;;) {
        chunk_t *slot = ring_get_putslot(io->ring);
#if defined(ZLIB)
        if (gzf) {
            int n = gzread(gzf, slot->buf, (unsigned)slot->cap);
            if (n <= 0) { slot->sz = 0; ring_commit_put(io->ring); break; }
            slot->sz = (size_t)n;
            ring_commit_put(io->ring);
            continue;
        }
#endif
        size_t n = fread(slot->buf, 1, slot->cap, fp);
        if (n == 0) { slot->sz = 0; ring_commit_put(io->ring); break; }
        slot->sz = n;
        ring_commit_put(io->ring);
    }

#if defined(ZLIB)
    if (gzf) gzclose(gzf);
#endif
    if (fp && fp!=stdin) fclose(fp);

    pthread_mutex_lock(&io->ring->mu);
    io->ring->closed = 1;
    pthread_cond_broadcast(&io->ring->cv_get);
    pthread_mutex_unlock(&io->ring->mu);
    return NULL;
}

// --- Public API ---------------------------------------------------------------
int kssd3_sketch_fast(const kssd3_sketch_cfg *cfg, kssd3_emit_fn emit, void *emit_user) {
    if (!cfg || !emit) return 1;
    if (cfg->k <= 0 || cfg->k > 31) { fprintf(stderr, "k must be 1..31\n"); return 2; }
    if (cfg->c < 2)                  { fprintf(stderr, "c must be >=2\n"); return 2; }
    int threads = cfg->threads > 0 ? cfg->threads : 1;
    size_t shard_bits = cfg->shard_bits ? cfg->shard_bits : 8; // 256 shards
    size_t ring_chunks = (size_t)threads * 2 + 2;
    size_t io_mb = cfg->io_buf_mb ? cfg->io_buf_mb : 4;

    // Init keep rule
    fmh_t fmh; fmh_init(&fmh, cfg->c);

    // Ring
    ring_t ring; ring_init(&ring, ring_chunks, io_mb);

    // I/O thread
    io_ctx_t io = { .ring = &ring, .cfg = cfg, .err = 0 };
    pthread_t io_thr;
    pthread_create(&io_thr, NULL, io_main, &io);

    // Workers
    worker_ctx_t **workers = (worker_ctx_t**)calloc((size_t)threads, sizeof(worker_ctx_t*));
    pthread_t *tids = (pthread_t*)calloc((size_t)threads, sizeof(pthread_t));
    if (!workers || !tids) DIE("OOM workers");

    for (int t=0;t<threads;t++) {
        workers[t] = (worker_ctx_t*)calloc(1, sizeof(worker_ctx_t));
        workers[t]->cfg = cfg;
        workers[t]->fmh = fmh;
        workers[t]->ring = &ring;
        workers[t]->shard_bits = shard_bits;
        pthread_create(&tids[t], NULL, worker_main, workers[t]);
    }

    // Join I/O
    pthread_join(io_thr, NULL);

    // Join workers
    for (int t=0;t<threads;t++) pthread_join(tids[t], NULL);

    // Emit
    emit_ctx_t ectx = { .nthreads = (size_t)threads, .workers = workers, .emit = emit, .emit_user = emit_user };
    emit_all(&ectx);

    // Cleanup
    for (int t=0;t<threads;t++) {
        shardset_t *ss = &workers[t]->ss;
        for (size_t s=0;s<ss->shards;s++) map_free(&ss->maps[s]);
        free(ss->maps);
        free(workers[t]);
    }
    free(workers);
    free(tids);
    ring_free(&ring);
    return 0;
}

// TSV helper
static void tsv_emit(uint64_t key, uint32_t count, void *user) {
    FILE *fp = (FILE*)user;
    fprintf(fp, "%" PRIu64 "\t%u\n", key, count);
}

int kssd3_sketch_fast_to_tsv(const kssd3_sketch_cfg *cfg, void *file_ptr) {
    FILE *fp = (FILE*)file_ptr;
    if (!fp) fp = stdout;
    return kssd3_sketch_fast(cfg, tsv_emit, fp);
}
