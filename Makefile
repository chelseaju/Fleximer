## Makefile

CC=g++

LIBINCLUDE=third-party/include/
LIB=third-party/lib/


CFLAGS=-I$(LIBINCLUDE) -std=c++11
LDFLAGS=-L$(LIB) -lsdsl -ldivsufsort -ldivsufsort64 --static

all: clean sigmer_generation sigmer_selection sigmer_count EM

sigmer_generation:
	$(CC) $(CFLAGS) sigmer_generation.cpp HELPER.cpp SEQCLUSTER.cpp -Ofast -march=native -o sigmer_generation $(LDFLAGS)

sigmer_selection:
	$(CC) $(CFLAGS) sigmer_selection.cpp HELPER.cpp SEQCLUSTER.cpp -Ofast -march=native -o sigmer_selection $(LDFLAGS)

sigmer_count:
	$(CC) $(CFLAGS) sigmer_count.cpp HELPER.cpp SEQCLUSTER.cpp aho01new.cpp -Ofast -march=native -o sigmer_count $(LDFLAGS)

EM:
	$(CC) $(CFLAGS) EM.cpp HELPER.cpp SEQCLUSTER.cpp -Ofast -march=native -o EM $(LDFLAGS)


clean:
	rm -rf *.o  sigmer_generation sigmer_selection sigmer_count EM
