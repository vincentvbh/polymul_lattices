
#include <stdio.h>
#include <stdint.h>

#define Q (8380417)

int32_t reduce_Dilithium(int32_t a){
    return a - ( (a + (1U << 22)) >> 23) * Q;
}

#define BOUNDMIN (-4186113)
#define BOUNDMAX (4194303)

int main(void){

    int32_t min, max;

    int32_t a;

    min = BOUNDMIN;
    max = BOUNDMAX;

    for(int32_t i = 0; i < Q; i++){
        a = reduce_Dilithium(i);
        if(max < a){
            max = a;
        }
        if(min > a){
            min = a;
        }
    }

    printf("%d\n%d\n", min, max);

}

