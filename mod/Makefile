
CC ?= gcc

CFLAGS += -mcpu=native -mtune=native -march=native -O3

COMMON_PATH = .

CFLAGS += -I$(COMMON_PATH)

COMMON_SOURCE =

ASM_HEADERs =
ASM_SOURCEs =
C_HEADERs =
C_SOURCEs =

SOURCEs = $(ASM_SOURCEs) $(C_SOURCEs) $(COMMON_SOURCE)
HEADERs = $(ASM_HEADERs) $(C_HEADERs)

all: dilithium_range

dilithium_range: dilithium_range.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@


.PHONY: clean
clean:
	rm -f dilithium_range



