CC = mpicc
CCFLAGS = -g -std=c99 -O3
LIBS = -lm

BINS = mc

all: $(BINS)


%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)
