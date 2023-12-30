
#include <stdint.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <assert.h>

#include "tools.h"
#include "naive_mult.h"
#include "gen_table.h"
#include "ntt_c.h"

#define ARRAY_N 1536

#define Q (7681)

int16_t mod = Q;

void memberZ(void *des, void *src){
    cmod_int16(des, src, &mod);
}

void addZ(void *des, void *src1, void *src2){
    addmod_int16(des, src1, src2, &mod);
}

void subZ(void *des, void *src1, void *src2){
    submod_int16(des, src1, src2, &mod);
}

void mulZ(void *des, void *src1, void *src2){
    mulmod_int16(des, src1, src2, &mod);
}

void expZ(void *des, void *src, size_t e){
    expmod_int16(des, src, e, &mod);
}

struct commutative_ring coeff_ring = {
    .sizeZ = sizeof(int16_t),
    .memberZ = memberZ,
    .addZ = addZ,
    .subZ = subZ,
    .mulZ = mulZ,
    .expZ = expZ
};

void memberZ_convol(void *des, void *src){

    size_t size = sizeof(int16_t);

    for(size_t i = 0; i < 3; i++){
        cmod_int16(des + i * size, src + i * size, &mod);
    }

}

void addZ_convol(void *des, void *src1, void *src2){

    size_t size = sizeof(int16_t);

    for(size_t i = 0; i < 3; i++){
        addmod_int16(des + i * size, src1 + i * size, src2 + i * size, &mod);
    }

}

void subZ_convol(void *des, void *src1, void *src2){

    size_t size = sizeof(int16_t);

    for(size_t i = 0; i < 3; i++){
        submod_int16(des + i * size, src1 + i * size, src2 + i * size, &mod);
    }

}

void mulZ_convol(void *des, void *src1, void *src2){

    size_t twiddle = 1;

    naive_mulR(des,
        src1, src2, 3, &twiddle, coeff_ring);

}

void expZ_convol(void *des, void *src, size_t e){
    return;
}

struct commutative_ring convol_ring = {
    .sizeZ = sizeof(int16_t) * 3,
    .memberZ = memberZ_convol,
    .addZ = addZ_convol,
    .subZ = subZ_convol,
    .mulZ = mulZ_convol,
    .expZ = expZ_convol
};

struct compress_profile profile;

#define BUFF_MAX 4096

int16_t bufflo[BUFF_MAX];
int16_t buffhi[BUFF_MAX];

int main(void){

    int16_t poly1[ARRAY_N], poly2[ARRAY_N];
    int16_t ref[ARRAY_N], res[ARRAY_N];

    int16_t poly1_NTT[ARRAY_N], poly2_NTT[ARRAY_N];
    int16_t res_NTT[ARRAY_N];

    int16_t twiddle_convol[3];

    int16_t omega, zeta, twiddle, scale, t;

    for(size_t i = 0; i < ARRAY_N; i++){
        t = rand();
        cmod_int16(poly1 + i, &t, &mod);
        t = rand();
        cmod_int16(poly2 + i, &t, &mod);
    }

    twiddle = 1;
    naive_mulR(ref,
        poly1, poly2, ARRAY_N, &twiddle, coeff_ring);

    for(size_t i = 0; i < ARRAY_N; i++){
        poly1_NTT[(i % 512) * 3 + (i % 3)] = poly1[i];
        poly2_NTT[(i % 512) * 3 + (i % 3)] = poly2[i];
    }

    twiddle_convol[0] = 1;
    twiddle_convol[1] = 0;
    twiddle_convol[2] = 0;
    naive_mulR(res_NTT,
        poly1_NTT, poly2_NTT, 512, &twiddle_convol, convol_ring);

    for(size_t i = 0; i < ARRAY_N; i++){
        res[i] = res_NTT[(i % 512) * 3 + (i % 3)];
    }

    for(size_t i = 0; i < ARRAY_N; i++){
        // if(ref[i] != res[i]){
        //     printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        // }
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");

}











