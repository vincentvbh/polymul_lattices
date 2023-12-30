
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

#define ARRAY_N 64
#define NTT_N 64
#define LOGNTT_N 6

#define Q (65537)

// #define OMEGA (2)
// #define OMEGA_INV (-32768)
#define OMEGA (-4080)
#define OMEGA_INV (-2040)

int32_t mod = Q;

void memberZ(void *des, void *src){
    cmod_int32(des, src, &mod);
}

void addZ(void *des, void *src1, void *src2){
    addmod_int32(des, src1, src2, &mod);
}

void subZ(void *des, void *src1, void *src2){
    submod_int32(des, src1, src2, &mod);
}

void mulZ(void *des, void *src1, void *src2){
    mulmod_int32(des, src1, src2, &mod);
}

void expZ(void *des, void *src, size_t e){
    expmod_int32(des, src, e, &mod);
}

struct commutative_ring coeff_ring = {
    .sizeZ = sizeof(int32_t),
    .memberZ = memberZ,
    .addZ = addZ,
    .subZ = subZ,
    .mulZ = mulZ,
    .expZ = expZ
};

struct compress_profile profile;

#define BUFF_MAX 4096

int32_t bufflo[BUFF_MAX];
int32_t buffhi[BUFF_MAX];

int32_t __sq[1u << 17];
int32_t __sqrt[1u << 17];

int main(void){

    int32_t poly1[ARRAY_N], poly2[ARRAY_N];
    int32_t ref[ARRAY_N], res[ARRAY_N];

    int32_t omega, zeta, twiddle, scale, t;

    for(size_t i = 0; i < ARRAY_N; i++){
        t = rand();
        cmod_int32(poly1 + i, &t, &mod);
        t = rand();
        cmod_int32(poly2 + i, &t, &mod);
    }

    twiddle = 1;
    naive_mulR(ref,
        poly1, poly2, ARRAY_N, &twiddle, coeff_ring);

    for(size_t i = 0; i < Q; i++){
        __sq[i] = __sqrt[i] = -1;
    }

    for(size_t i = 0; i < Q; i++){
        t = i;
        mulmod_int32(__sq + i, &t, &t, &mod);
    }

    for(size_t i = 0; i < Q; i++){
        t = __sq[i];
        if(t < 0){
            t += Q;
        }
        __sqrt[t] = i;
        cmod_int32(__sqrt + t, __sqrt + t, &mod);
    }

    printf("%d\n", __sqrt[2]);


    profile.array_n = ARRAY_N;
    profile.ntt_n = NTT_N;
    profile.log_ntt_n = LOGNTT_N;

    profile.compressed_layers = LOGNTT_N;
    for(size_t i = 0; i < profile.compressed_layers; i++){
        profile.merged_layers[i] = 1;
    }

    zeta = 1;
    omega = OMEGA;
    scale = 1;
    gen_streamlined_DWT_table(buffhi,
        &scale, &omega, &zeta, profile, 0, coeff_ring);

    for(size_t i = 0; i < 63; i++){
        printf("%d, ", buffhi[i]);
    }
    printf("\n");

    compressed_CT_NTT(poly1,
        0, LOGNTT_N - 1, buffhi, profile, coeff_ring);
    compressed_CT_NTT(poly2,
        0, LOGNTT_N - 1, buffhi, profile, coeff_ring);

    point_mul(res, poly1, poly2, ARRAY_N, 1, coeff_ring);

    zeta = 1;
    omega = OMEGA_INV;
    scale = 1;
    gen_streamlined_DWT_table(buffhi,
        &scale, &omega, &zeta, profile, 0, coeff_ring);

    compressed_GS_NTT(res,
        0, LOGNTT_N - 1, buffhi, profile, coeff_ring);

    scale = NTT_N;
    for(size_t i = 0; i < ARRAY_N; i++){
        coeff_ring.mulZ(ref + i, ref + i, &scale);
    }

    for(size_t i = 0; i < ARRAY_N; i++){
        // if(ref[i] != res[i]){
        //     printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        // }
        assert(ref[i] == res[i]);
    }

    printf("Test finished!\n");

}








