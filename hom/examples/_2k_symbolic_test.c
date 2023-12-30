
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

#include "tools.h"
#include "naive_mult.h"

// x^256 - 1
// x^8 - y, y^32 - 1
// x^16 + 1, y^32 - 1
#define ARRAY_N 256
#define INNER_N 16

#define Q (1 << 13)

int32_t mod = Q;

#define COEFF_TYPE int32_t
#define COEFF_SIZE (sizeof(COEFF_TYPE))


void memberZ(void *des, void *src){
    *(COEFF_TYPE*)des = *(COEFF_TYPE*)src;
}

void addZ(void *des, void *src1, void *src2){
    *(COEFF_TYPE*)des = (*(COEFF_TYPE*)src1) + (*(COEFF_TYPE*)src2);
}

void subZ(void *des, void *src1, void *src2){
    *(COEFF_TYPE*)des = (*(COEFF_TYPE*)src1) - (*(COEFF_TYPE*)src2);
}

void mulZ(void *des, void *src1, void *src2){
    *(COEFF_TYPE*)des = (*(COEFF_TYPE*)src1) * (*(COEFF_TYPE*)src2);
}

void expZ(void *des, void *src, size_t e){
    return;
}

struct commutative_ring coeff_ring = {
    .sizeZ = COEFF_SIZE,
    .memberZ = memberZ,
    .addZ = addZ,
    .subZ = subZ,
    .mulZ = mulZ,
    .expZ = expZ
};

void memberZ_negacyclic(void *des, void *src){

    for(size_t i = 0; i < INNER_N; i++){
        coeff_ring.memberZ(des + i * COEFF_SIZE, src + i * COEFF_SIZE);
    }

}

void addZ_negacyclic(void *des, void *src1, void *src2){

    for(size_t i = 0; i < INNER_N; i++){
        coeff_ring.addZ(des + i * COEFF_SIZE, src1 + i * COEFF_SIZE, src2 + i * COEFF_SIZE);
    }

}

void subZ_negacyclic(void *des, void *src1, void *src2){

    for(size_t i = 0; i < INNER_N; i++){
        coeff_ring.subZ(des + i * COEFF_SIZE, src1 + i * COEFF_SIZE, src2 + i * COEFF_SIZE);
    }

}

void mulZ_negacyclic(void *des, void *src1, void *src2){

    COEFF_TYPE twiddle = -1;

    naive_mulR(des, src1, src2, INNER_N, &twiddle, coeff_ring);

}

void expZ_negacyclic(void *des, void *src, size_t e){
    return;
}

struct commutative_ring negacyclic_ring = {
    .sizeZ = COEFF_SIZE * INNER_N,
    .memberZ = memberZ_negacyclic,
    .addZ = addZ_negacyclic,
    .subZ = subZ_negacyclic,
    .mulZ = mulZ_negacyclic,
    .expZ = expZ_negacyclic
};

// 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
__attribute((aligned(32))) COEFF_TYPE twiddle_symbolic_CT_tabe[16][INNER_N] = {
{1, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 1, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 1, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 1, 0, 0, 0},

{0, 0, 1, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 1, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 1, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 1, 0},

{0, 1, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 1, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 1, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 1, 0, 0},

{0, 0, 0, 1,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 1,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 1,
 0, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 1}
};

// 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
__attribute((aligned(32))) COEFF_TYPE twiddle_symbolic_GS_tabe[16][INNER_N] = {
{1, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 -1, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 -1, 0, 0, 0},
{0, 0, 0, 0,
 -1, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},

{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, -1, 0},
{0, 0, 0, 0,
 0, 0, -1, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, -1, 0,
 0, 0, 0, 0},
{0, 0, -1, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},

{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, -1},
{0, 0, 0, 0,
 0, 0, 0, -1,
 0, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, -1,
 0, 0, 0, 0},
{0, 0, 0, -1,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},

{0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, -1, 0, 0},
{0, 0, 0, 0,
 0, -1, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0},
{0, 0, 0, 0,
 0, 0, 0, 0,
 0, -1, 0, 0,
 0, 0, 0, 0},
{0, -1, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0, 0, 0, 0}
};

void twiddle_symbolic_test(){

    COEFF_TYPE twiddle[INNER_N];
    bool correct;

    for(size_t i = 0; i < 16; i++){
        negacyclic_ring.mulZ(twiddle, twiddle_symbolic_CT_tabe[i], twiddle_symbolic_GS_tabe[i]);
        correct = 1;
        if(twiddle[0] != 1){
            correct = 0;
        }
        for(size_t j = 1; j < INNER_N; j++){
            if(twiddle[j] != 0){
                correct = 0;
            }
        }
        if(correct == 0){
            printf("%4zu is wrong!\n", i);
            for(size_t j = 0; j < INNER_N; j++){
                printf("%4zu: %12d\n", j, twiddle[j]);
            }
        }
    }

}



void FFT_symbolic(COEFF_TYPE *src){

    COEFF_TYPE buff[INNER_N];
    COEFF_TYPE twiddle[INNER_N];
    size_t jump;


    jump = 32;
    for(size_t i = 0; i < 1; i++){
        memmove(twiddle, twiddle_symbolic_CT_tabe[i], negacyclic_ring.sizeZ);
        for(size_t j = 0; j < (jump >> 1); j++){
            negacyclic_ring.mulZ(buff, src + (j + (jump >> 1)) * INNER_N, twiddle);
            negacyclic_ring.subZ(src + (j +  (jump >> 1)) * INNER_N, src + (j + 0) * INNER_N, buff);
            negacyclic_ring.addZ(src + (j +            0) * INNER_N, src + (j + 0) * INNER_N, buff);
        }
        src += jump * INNER_N;
    }
    src -= 32 * INNER_N;

    jump = 16;
    for(size_t i = 0; i < 2; i++){
        memmove(twiddle, twiddle_symbolic_CT_tabe[i], negacyclic_ring.sizeZ);
        for(size_t j = 0; j < (jump >> 1); j++){
            negacyclic_ring.mulZ(buff, src + (j + (jump >> 1)) * INNER_N, twiddle);
            negacyclic_ring.subZ(src + (j +  (jump >> 1)) * INNER_N, src + (j + 0) * INNER_N, buff);
            negacyclic_ring.addZ(src + (j +            0) * INNER_N, src + (j + 0) * INNER_N, buff);
        }
        src += jump * INNER_N;
    }
    src -= 32 * INNER_N;

    jump = 8;
    for(size_t i = 0; i < 4; i++){
        memmove(twiddle, twiddle_symbolic_CT_tabe[i], negacyclic_ring.sizeZ);
        for(size_t j = 0; j < (jump >> 1); j++){
            negacyclic_ring.mulZ(buff, src + (j + (jump >> 1)) * INNER_N, twiddle);
            negacyclic_ring.subZ(src + (j +  (jump >> 1)) * INNER_N, src + (j + 0) * INNER_N, buff);
            negacyclic_ring.addZ(src + (j +            0) * INNER_N, src + (j + 0) * INNER_N, buff);
        }
        src += jump * INNER_N;
    }
    src -= 32 * INNER_N;

    jump = 4;
    for(size_t i = 0; i < 8; i++){
        memmove(twiddle, twiddle_symbolic_CT_tabe[i], negacyclic_ring.sizeZ);
        for(size_t j = 0; j < (jump >> 1); j++){
            negacyclic_ring.mulZ(buff, src + (j + (jump >> 1)) * INNER_N, twiddle);
            negacyclic_ring.subZ(src + (j +  (jump >> 1)) * INNER_N, src + (j + 0) * INNER_N, buff);
            negacyclic_ring.addZ(src + (j +            0) * INNER_N, src + (j + 0) * INNER_N, buff);
        }
        src += jump * INNER_N;
    }
    src -= 32 * INNER_N;

    jump = 2;
    for(size_t i = 0; i < 16; i++){
        memmove(twiddle, twiddle_symbolic_CT_tabe[i], negacyclic_ring.sizeZ);
        for(size_t j = 0; j < (jump >> 1); j++){
            negacyclic_ring.mulZ(buff, src + (j + (jump >> 1)) * INNER_N, twiddle);
            negacyclic_ring.subZ(src + (j +  (jump >> 1)) * INNER_N, src + (j + 0) * INNER_N, buff);
            negacyclic_ring.addZ(src + (j +            0) * INNER_N, src + (j + 0) * INNER_N, buff);
        }
        src += jump * INNER_N;
    }
    src -= 32 * INNER_N;

}

void iFFT_symbolic(COEFF_TYPE *src){

    COEFF_TYPE buff[INNER_N];
    COEFF_TYPE twiddle[INNER_N];
    size_t jump;

    jump = 2;
    for(size_t i = 0; i < 16; i++){
        memmove(twiddle, twiddle_symbolic_GS_tabe[i], negacyclic_ring.sizeZ);
        for(size_t j = 0; j < (jump >> 1); j++){
            negacyclic_ring.subZ(                    buff, src + (j +  0) * INNER_N, src + (j + (jump >> 1)) * INNER_N);
            negacyclic_ring.addZ(src + (j +  0) * INNER_N, src + (j +  0) * INNER_N, src + (j + (jump >> 1)) * INNER_N);
            negacyclic_ring.mulZ(src + (j + (jump >> 1)) * INNER_N, buff, twiddle);
        }
        src += jump * INNER_N;
    }
    src -= 32 * INNER_N;

    jump = 4;
    for(size_t i = 0; i < 8; i++){
        memmove(twiddle, twiddle_symbolic_GS_tabe[i], negacyclic_ring.sizeZ);
        for(size_t j = 0; j < (jump >> 1); j++){
            negacyclic_ring.subZ(                    buff, src + (j +  0) * INNER_N, src + (j + (jump >> 1)) * INNER_N);
            negacyclic_ring.addZ(src + (j +  0) * INNER_N, src + (j +  0) * INNER_N, src + (j + (jump >> 1)) * INNER_N);
            negacyclic_ring.mulZ(src + (j + (jump >> 1)) * INNER_N, buff, twiddle);
        }
        src += jump * INNER_N;
    }
    src -= 32 * INNER_N;

    jump = 8;
    for(size_t i = 0; i < 4; i++){
        memmove(twiddle, twiddle_symbolic_GS_tabe[i], negacyclic_ring.sizeZ);
        for(size_t j = 0; j < (jump >> 1); j++){
            negacyclic_ring.subZ(                    buff, src + (j +  0) * INNER_N, src + (j + (jump >> 1)) * INNER_N);
            negacyclic_ring.addZ(src + (j +  0) * INNER_N, src + (j +  0) * INNER_N, src + (j + (jump >> 1)) * INNER_N);
            negacyclic_ring.mulZ(src + (j + (jump >> 1)) * INNER_N, buff, twiddle);
        }
        src += jump * INNER_N;
    }
    src -= 32 * INNER_N;

    jump = 16;
    for(size_t i = 0; i < 2; i++){
        memmove(twiddle, twiddle_symbolic_GS_tabe[i], negacyclic_ring.sizeZ);
        for(size_t j = 0; j < (jump >> 1); j++){
            negacyclic_ring.subZ(                    buff, src + (j +  0) * INNER_N, src + (j + (jump >> 1)) * INNER_N);
            negacyclic_ring.addZ(src + (j +  0) * INNER_N, src + (j +  0) * INNER_N, src + (j + (jump >> 1)) * INNER_N);
            negacyclic_ring.mulZ(src + (j + (jump >> 1)) * INNER_N, buff, twiddle);
        }
        src += jump * INNER_N;
    }
    src -= 32 * INNER_N;

    jump = 32;
    for(size_t i = 0; i < 1; i++){
        memmove(twiddle, twiddle_symbolic_GS_tabe[i], negacyclic_ring.sizeZ);
        for(size_t j = 0; j < (jump >> 1); j++){
            negacyclic_ring.subZ(                    buff, src + (j +  0) * INNER_N, src + (j + (jump >> 1)) * INNER_N);
            negacyclic_ring.addZ(src + (j +  0) * INNER_N, src + (j +  0) * INNER_N, src + (j + (jump >> 1)) * INNER_N);
            negacyclic_ring.mulZ(src + (j + (jump >> 1)) * INNER_N, buff, twiddle);
        }
        src += jump * INNER_N;
    }
    src -= 32 * INNER_N;

}

int main(void){

    int32_t poly1[ARRAY_N], poly2[ARRAY_N];
    int32_t poly1_NTT[32 * INNER_N], poly2_NTT[32 * INNER_N];
    int32_t res_NTT[32 * INNER_N];
    int32_t ref[ARRAY_N], res[ARRAY_N];

    int32_t omega, twiddle, scale, t;

    int32_t twiddle_negacyclic[INNER_N];
    int32_t buff[INNER_N];

    for(size_t i = 0; i < ARRAY_N; i++){
        t = rand();
        cmod_int32(poly1 + i, &t, &mod);
        t = rand();
        cmod_int32(poly2 + i, &t, &mod);
    }

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, ARRAY_N, &twiddle, coeff_ring);
    for(size_t i = 0; i < ARRAY_N; i++){
        cmod_int32(ref + i, ref + i, &mod);
    }

    for(size_t i = 0; i < 32; i++){
        for(size_t j = 0; j < 8; j++){
            poly1_NTT[i * 16 + j] = poly1[i * 8 + j];
            poly2_NTT[i * 16 + j] = poly2[i * 8 + j];
        }
        for(size_t j = 8; j < 16; j++){
            poly1_NTT[i * 16 + j] = 0;
            poly2_NTT[i * 16 + j] = 0;
        }
    }

    for(size_t i = 0; i < 32; i++){
        for(size_t j = 0; j < 16; j++){
            res_NTT[i * 16 + j] = 0;
        }
    }

    twiddle_negacyclic[0] = 1;
    for(size_t i = 1; i < 16; i++){
        twiddle_negacyclic[i] = 0;
    }

#if 1

    twiddle_symbolic_test();

    FFT_symbolic(poly1_NTT);
    FFT_symbolic(poly2_NTT);

    for(size_t i = 0; i < 32; i++){
        negacyclic_ring.mulZ(res_NTT + i * INNER_N, poly1_NTT + i * INNER_N, poly2_NTT + i * INNER_N);
    }

    iFFT_symbolic(res_NTT);

    for(size_t i = 0; i < 32 * INNER_N; i++){
        res_NTT[i] >>= 5;
    }

#else

    naive_mulR(res_NTT, poly1_NTT, poly2_NTT, 32, &twiddle_negacyclic, negacyclic_ring);

#endif

    for(size_t i = 1; i < 32; i++){
        for(size_t j = 0; j < 8; j++){
            coeff_ring.addZ(res_NTT + i * 16 + j, res_NTT + i * 16 + j, res_NTT + (i - 1) * 16 + j + 8);
        }
    }

    for(size_t j = 0; j < 8; j++){
        coeff_ring.addZ(res_NTT + 0 * 16 + j, res_NTT + 0 * 16 + j, res_NTT + 31 * 16 + j + 8);
    }

    for(size_t i = 0; i < 32; i++){
        for(size_t j = 0; j < 8; j++){
            res[i * 8 + j] = res_NTT[i * 16 + j];
        }
    }

    for(size_t i = 0; i < ARRAY_N; i++){
        cmod_int32(res + i, res + i, &mod);
    }

    for(size_t i = 0; i < ARRAY_N; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12x, %12x\n", i, ref[i], res[i]);
        }
    }

    printf("test finished!\n");

}









