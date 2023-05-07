CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

BIN = matmul

all: $(BIN)

matmul: test_arrays_revised.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
