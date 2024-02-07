
CC = gcc

CFLAGS += -O0 -Wall

# $(warning $(CFLAGS))

COMMON_PATH = ../common

CFLAGS += -I$(COMMON_PATH)

GMP_PATH = /opt/homebrew/opt/gmp

GMP_FLAGS = -I$(GMP_PATH)/include -L$(GMP_PATH)/lib -lgmp

COMMON_SOURCE = $(COMMON_PATH)/tools.c $(COMMON_PATH)/naive_mult.c $(COMMON_PATH)/gen_table.c $(COMMON_PATH)/ntt_c.c

ASM_HEADERs =
ASM_SOURCEs =
C_HEADERs =
C_SOURCEs =

SOURCEs = $(ASM_SOURCEs) $(C_SOURCEs) $(COMMON_SOURCE)
HEADERs = $(ASM_HEADERs) $(C_HEADERs)

all: DWT FNT GT Nussbaumer Schoenhage TMVP Karatsuba BigIntMul

DWT: DWT.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

FNT: FNT.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

GT: GT.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Nussbaumer: Nussbaumer.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Schoenhage: Schoenhage.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

TMVP: TMVP.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Karatsuba: Karatsuba.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

BigIntMul: BigIntMul.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $(GMP_FLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f DWT
	rm -f FNT
	rm -f GT
	rm -f Nussbaumer
	rm -f Schoenhage
	rm -f TMVP
	rm -f Karatsuba
	rm -f BigIntMul



