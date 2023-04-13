CC = mpicc
CCFLAGS = -g -O3 -Wall
LIBS = -lm

BINS = main

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

