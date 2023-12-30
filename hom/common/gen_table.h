
#ifndef GEN_TABLE_H
#define GEN_TABLE_H

#include <stddef.h>
#include <stdbool.h>

#include "tools.h"

// addZ, mulZ, expZ

// ================================
// generate twiddle factors for cyclic NTT with Cooley-Tukey butterflies

void gen_CT_table(
    void *des,
    void *scale, void *omega,
    struct compress_profile _profile,
    struct commutative_ring ring
    );

// ================================
// generate twiddle factors for DWT with Cooley-Tukey butterflies

void gen_DWT_table(
    void *des,
    void *scale, void *omega, void *zeta,
    struct compress_profile _profile,
    struct commutative_ring ring
    );

// ================================
// generate twiddle factors for cyclic iNTT with Cooley-Tukey butterflies

void gen_inv_CT_table(
    void *des,
    void *scale, void *omega,
    struct compress_profile _profile,
    struct commutative_ring ring
    );

// ================================
// generate twiddle factors for DWT with Cooley-Tukey butterflies
// the table is re-ordered according to _profile

void gen_streamlined_DWT_table(
    void *des,
    void *scale, void *omega, void *zeta,
    struct compress_profile _profile, bool pad,
    struct commutative_ring ring
    );

// ================================
// generate twiddle factors for cyclic iNTT with Cooley-Tukey butterflies
// the table is re-ordered according to _profile

void gen_streamlined_inv_CT_table(
    void *des,
    void *scale, void *omega,
    struct compress_profile _profile, bool pad,
    struct commutative_ring ring
    );

// ================================
// generate twiddle factors for twisting (x^NTT_N - omega^NTT_N) to (x^NTT_N - 1)

void gen_twist_table(
    void *des,
    void *scale, void *omega,
    struct compress_profile _profile,
    struct commutative_ring ring
    );

// ================================
// generate twiddle factors for multiplication in
// x^(ARRAY_N / NTT_N) +- omega^i

void gen_mul_table(
    void *des,
    void *scale, void *omega,
    struct compress_profile _profile,
    struct commutative_ring ring
    );

#endif

