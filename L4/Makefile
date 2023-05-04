CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

BIN = E1

all: $(BIN)

E1: E1.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
