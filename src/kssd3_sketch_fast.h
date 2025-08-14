// kssd3_sketch_fast.h
// High-throughput FracMinHash sketcher core for FASTQ/FASTA (plain or .gz).
// No CLI; meant to be called from your existing command_sketch.c.

#ifndef KSSD3_SKETCH_FAST_H
#define KSSD3_SKETCH_FAST_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Emit callback: called once per (key,count) at the end (after per-thread merges)
typedef void (*kssd3_emit_fn)(uint64_t key, uint32_t count, void *user);

// Configuration for the sketcher.
typedef struct {
    const char *input_path;     // "-" for stdin; supports .gz if built with zlib/libdeflate
    int         k;              // k-mer length (supports up to 31 efficiently in 64 bits)
    uint32_t    c;              // FracMinHash denominator; if power-of-two, keep test is bitmask
    int         threads;        // worker threads (>=1). One extra I/O thread is created internally.
    int         is_fastq;       // 1 = FASTQ, 0 = FASTA (best-effort autodetect if <0)
    size_t      shard_bits;     // number of top hash bits for sharding (e.g., 8 => 256 shards). 8–10 is good.
    size_t      io_buf_mb;      // size of I/O buffer in MB (e.g., 4–16MB)
} kssd3_sketch_cfg;

// Main API: run fast sketch and emit (key,count) via callback.
// Returns 0 on success; nonzero on error.
int kssd3_sketch_fast(const kssd3_sketch_cfg *cfg, kssd3_emit_fn emit, void *emit_user);

// Convenience helper: write results as TSV to a FILE*.
// Columns: hash\tcount\n
int kssd3_sketch_fast_to_tsv(const kssd3_sketch_cfg *cfg, void *file_ptr /* FILE* */);

#ifdef __cplusplus
}
#endif

#endif // KSSD3_SKETCH_FAST_H
