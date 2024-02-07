

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "tools.h"
#include "naive_mult.h"

#include "gmp.h"

// ================
// We compute the product of two size-128 polynomials in Z_{2^64}[x] via
// integer-to-polynomial reduction.

// ================
// Optimization guide.
/*

1. Figure out what is the fastest approach multiplying size-128 polynomials in Z_Q[x]
   with Q >= (2^16)^2 * (2048 / 16) = 2^39. We compute entirely over Z_{2^64}. Figure out
   if there are better choices.

*/

#define LIMBS 32
#define ARRAY_N 128

// ================
// Z_{2^64}

void memberZ(void *des, void *src){
    *(uint64_t*)des = *(uint64_t*)src;
}

void addZ(void *des, void *src1, void *src2){
    *(uint64_t*)des = (*(uint64_t*)src1) + (*(uint64_t*)src2);
}

void subZ(void *des, void *src1, void *src2){
    *(uint64_t*)des = (*(uint64_t*)src1) - (*(uint64_t*)src2);
}

void mulZ(void *des, void *src1, void *src2){
    *(uint64_t*)des = (*(uint64_t*)src1) * (*(uint64_t*)src2);
}

void expZ(void *des, void *src, size_t e){

    uint64_t src_v = *(uint64_t*)src;
    uint64_t tmp_v;

    tmp_v = 1;
    for(; e; e >>= 1){
        if(e & 1){
            tmp_v = tmp_v * src_v;
        }
        src_v = src_v * src_v;
    }

    memmove(des, &tmp_v, sizeof(uint64_t));
}

struct commutative_ring coeff_ring = {
    .sizeZ = sizeof(uint64_t),
    .memberZ = memberZ,
    .addZ = addZ,
    .subZ = subZ,
    .mulZ = mulZ,
    .expZ = expZ
};

void gmp_mul(uint64_t *res, const uint64_t *src1, const uint64_t *src2, size_t len){

    mpz_t x, y, z;

    mpz_init(x);
    mpz_init(y);
    mpz_init(z);

    mpz_import(x, len, -1, sizeof(uint64_t), -1, 0, src1);
    mpz_import(y, len, -1, sizeof(uint64_t), -1, 0, src2);

    mpz_mul(z, x, y);

    memset(res, 0, 2 * len * sizeof(uint64_t));
    mpz_export(res, NULL, -1, sizeof(uint64_t), -1, 0, z);

    mpz_clears(x, y, z, NULL);

}

void bits_chunk16(uint64_t *des, uint64_t *src, size_t len){

    uint16_t *ptr;

    ptr = (uint16_t*)src;
    for(size_t i = 0; i < 4 * len; i++){
        des[i] = (uint64_t)(*ptr++);
    }

}

void bits_dechunk16(uint64_t *des, uint64_t *src, size_t len){

    uint64_t a, t;
    uint16_t *ptr;

    ptr = (uint16_t*)des;
    a = 0;
    for(size_t i = 0; i < len; i++){
        t = *src++;
        t += a;
        *ptr++ = (uint16_t)t;
        a = t >> 16;
    }

    *ptr++ = (uint16_t)(a);

}

int main(void){

    uint64_t src1[LIMBS], src2[LIMBS];
    uint64_t ref[2 * LIMBS], res[2 * LIMBS];
    uint64_t poly1[ARRAY_N], poly2[ARRAY_N];
    uint64_t poly_prod[2 * ARRAY_N];

    for(size_t i = 0; i < LIMBS; i++){
        src1[i] = (((uint64_t)rand()) << 32 ) | rand();
        src2[i] = (((uint64_t)rand()) << 32 ) | rand();
    }

    gmp_mul(ref, src1, src2, LIMBS);

    printf("gmp_mul finished!\n");

    bits_chunk16(poly1, src1, LIMBS);
    bits_chunk16(poly2, src2, LIMBS);

    memset(poly_prod, 0, 2 * ARRAY_N * sizeof(uint64_t));
    naive_mul_long(poly_prod, poly1, poly2, ARRAY_N, coeff_ring);

    bits_dechunk16(res, poly_prod, 2 * ARRAY_N - 1);

    for(size_t i = 0; i < 2 * LIMBS; i++){
        assert(ref[i] == res[i]);
    }

    printf("naive_mul_long finished!\n");

}








