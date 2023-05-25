CC = mpicc
CFLAGS = -std=c99 -g -Wall -O3 -ffast-math
LIBS = -lm

BIN = cg

all: $(BIN)

cg: main.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
