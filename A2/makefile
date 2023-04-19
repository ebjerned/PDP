CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

BIN = matrix

all: $(BIN)

matrix: test_arrays.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
