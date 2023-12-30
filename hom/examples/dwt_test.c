
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

#define ARRAY_N 512
#define NTT_N 512
#define LOGNTT_N 9

#define Q (12289)

#define OMEGA (49)
#define OMEGA_INV (1254)

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

struct compress_profile profile;

#define BUFF_MAX 4096

int16_t bufflo[BUFF_MAX];
int16_t buffhi[BUFF_MAX];


int main(void){

    int16_t poly1[ARRAY_N], poly2[ARRAY_N];
    int16_t ref[ARRAY_N], res[ARRAY_N];

    int16_t omega, zeta, twiddle, scale, t;

    for(size_t i = 0; i < ARRAY_N; i++){
        t = rand();
        cmod_int16(poly1 + i, &t, &mod);
        t = rand();
        cmod_int16(poly2 + i, &t, &mod);
    }

    twiddle = -1;
    naive_mulR(ref,
        poly1, poly2, ARRAY_N, &twiddle, coeff_ring);

    profile.array_n = ARRAY_N;
    profile.ntt_n = NTT_N;
    profile.log_ntt_n = LOGNTT_N;

    profile.compressed_layers = 4;
    profile.merged_layers[0] = 3;
    profile.merged_layers[1] = 2;
    profile.merged_layers[2] = 2;
    profile.merged_layers[3] = 2;

    zeta = OMEGA;
    coeff_ring.expZ(&omega, &zeta, 2);
    scale = 1;
    gen_streamlined_DWT_table(buffhi,
        &scale, &omega, &zeta, profile, 0, coeff_ring);

    compressed_CT_NTT(poly1,
        0, 3, buffhi, profile, coeff_ring);
    compressed_CT_NTT(poly2,
        0, 3, buffhi, profile, coeff_ring);

    point_mul(res, poly1, poly2, ARRAY_N, 1, coeff_ring);

    zeta = OMEGA_INV;
    coeff_ring.expZ(&omega, &zeta, 2);
    scale = 1;
    gen_streamlined_DWT_table(buffhi,
        &scale, &omega, &zeta, profile, 0, coeff_ring);

    compressed_GS_NTT(res,
        0, 3, buffhi, profile, coeff_ring);

    scale = 512;
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








