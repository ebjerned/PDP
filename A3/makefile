CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

BIN = sort

all: $(BIN)

sort: sort.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)

