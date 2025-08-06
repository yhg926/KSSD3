#ifndef DNA_POPCOUNT_H
#define DNA_POPCOUNT_H

#include <stdint.h>

/* ---- Fallback popcount for 32-bit ---- */
static inline int popcount32_fallback(uint32_t x) {
    x = x - ((x >> 1) & 0x55555555U);
    x = (x & 0x33333333U) + ((x >> 2) & 0x33333333U);
    x = (x + (x >> 4)) & 0x0F0F0F0FU;
    x = x + (x >> 8);
    x = x + (x >> 16);
    return x & 0x3F;
}

/* ---- Fallback popcount for 64-bit ---- */
static inline int popcount64_fallback(uint64_t x) {
    x = x - ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    x = x + (x >> 8);
    x = x + (x >> 16);
    x = x + (x >> 32);
    return x & 0x7F;
}

/* ---- Portable popcount wrapper ---- */
static inline int popcount32(uint32_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcount(x);
#else
    return popcount32_fallback(x);
#endif
}

static inline int popcount64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    return popcount64_fallback(x);
#endif
}

/* ---- Compare two 2-bit encoded base sequences (uint32_t, ≤16 bases) ---- */
static inline int dna_popcount(uint32_t xXORy) {
    xXORy |= xXORy >> 1;
    xXORy &= 0x55555555U;
    return popcount32(xXORy);
}

/* ---- Compare two 2-bit encoded base sequences (uint64_t, ≤32 bases) ---- */
static inline int dna_popcount64(uint64_t xXORy) {
    xXORy |= xXORy >> 1;
    xXORy &= 0x5555555555555555ULL;
    return popcount64(xXORy);
}

#endif // BASE_DIFF_2BIT_H

