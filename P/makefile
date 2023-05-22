CC = mpicc
CFLAGS = -std=c99 -g -Wall -O3 -fopt-info-vec -fopt-info-vec-missed
LIBS = -lm

BIN = cg

all: $(BIN)

cg: test.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
