SHELL=/bin/bash
CC=mpicc
CFLAGS=-O2
LIBS=-lm
OUT_DIR=bin

.PHONY: directories all clean

all: directories mulRowCannon mulBlockColCannon

directories: $(OUT_DIR)

$(OUT_DIR):
	mkdir -p ./$@

mulRowCannon: ./src/mulRowCannon.c
	$(CC) $(CFLAGS) $< $(LIBS) -o ./bin/$@

mulBlockColCannon: ./src/mulBlockColCannon.c
	$(CC) $(CFLAGS) $< $(LIBS) -o ./bin/$@


clean:
	rm -rf $(OUT_DIR)
