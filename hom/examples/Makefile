
CC ?= gcc

CFLAGS += -mcpu=native -mtune=native -march=native -O3

COMMON_PATH = ../common

CFLAGS += -I$(COMMON_PATH)

COMMON_SOURCE = $(COMMON_PATH)/tools.c $(COMMON_PATH)/naive_mult.c $(COMMON_PATH)/gen_table.c $(COMMON_PATH)/ntt_c.c

ASM_HEADERs =
ASM_SOURCEs =
C_HEADERs =
C_SOURCEs =

SOURCEs = $(ASM_SOURCEs) $(C_SOURCEs) $(COMMON_SOURCE)
HEADERs = $(ASM_HEADERs) $(C_HEADERs)

all: dwt_test fnt_test gt_test symbolic_test _2k_symbolic_test

dwt_test: dwt_test.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

fnt_test: fnt_test.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

gt_test: gt_test.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

symbolic_test: symbolic_test.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

_2k_symbolic_test: _2k_symbolic_test.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@



.PHONY: clean
clean:
	rm -f dwt_test
	rm -f fnt_test
	rm -f gt_test
	rm -f symbolic_test
	rm -f _2k_symbolic_test



