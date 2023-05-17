CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm


BIN = quicksort

all: $(BIN)

quicksort: new_main.c	
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
clean:
	$(RM) $(BIN)

