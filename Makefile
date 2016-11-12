## Makefile

CC=g++
#LIBINCLUDE=/u/home/c/chelseaj/project/SDSL/include/
#LIBINCLUDE=/u/home/c/chelseaj/project/FunSkim/CSuffixSkim/third-party/include/
#LIBINCLUDE=/home/chelseaju/workspace/FunSkim/CSuffixSkim/third-party/include/
#LIB=/u/home/c/chelseaj/project/SDSL/lib/
#LIB=/u/home/c/chelseaj/project/FunSkim/CSuffixSkim/third-party/lib/
#LIB=/home/chelseaju/workspace/FunSkim/CSuffixSkim/third-party/lib

## QILIN
LIBINCLUDE=/home/chelseaju/FunSkim/FunSkim/CSuffixSkim/third-party/include/
LIB=/home/chelseaju/FunSkim/FunSkim/CSuffixSkim/third-party/lib/


CFLAGS=-I$(LIBINCLUDE) -std=c++11
LDFLAGS=-L$(LIB) -lsdsl -ldivsufsort -ldivsufsort64 --static

all: clean sigmer_generation sigmer_count_01_new EM

sigmer_generation:
	$(CC) $(CFLAGS) sigmer_generation_v10.cpp HELPER.cpp SEQCLUSTER.cpp -Ofast -march=native -o sigmer_generation $(LDFLAGS)

sigmer_selection:
	$(CC) $(CFLAGS) sigmer_selection_v92.cpp HELPER.cpp SEQCLUSTER.cpp -Ofast -march=native -o sigmer_selection $(LDFLAGS)

sigmer_cluster:
	$(CC) $(CFLAGS) sigmer_cluster.cpp HELPER.cpp SEQCLUSTER.cpp -o sigmer_cluster $(LDFLAGS)

sigmer_count_01_new:
	$(CC) $(CFLAGS) sigmer_count_v30.cpp HELPER.cpp SEQCLUSTER.cpp aho01new.cpp -Ofast -march=native -o sigmer_count_01_new $(LDFLAGS)

sigmer_count_01:
	$(CC) $(CFLAGS) sigmer_count_v22.cpp HELPER.cpp SEQCLUSTER.cpp aho01.cpp -Ofast -march=native -o sigmer_count_01 $(LDFLAGS)

sigmer_count_acgt:
	$(CC) $(CFLAGS) sigmer_count_v14.cpp HELPER.cpp SEQCLUSTER.cpp ahoacgt.cpp -Ofast -march=native -o sigmer_count_acgt $(LDFLAGS)

sigmer_estimation:
	$(CC) $(CFLAGS) sigmer_estimation_v6.cpp HELPER.cpp SEQCLUSTER.cpp -o sigmer_estimation $(LDFLAGS)

EM:
	$(CC) $(CFLAGS) EM_v15.cpp HELPER.cpp SEQCLUSTER.cpp -Ofast -march=native -o EM $(LDFLAGS)

trace_tree:
	$(CC) $(CFLAGS) trace_tree_v2.cpp HELPER.cpp SEQCLUSTER.cpp -o trace_tree $(LDFLAGS)

debug_sigmerIndex_to_binary:
	$(CC) $(CFLAGS) debug_sigmerIndex_to_binary.cpp HELPER.cpp SEQCLUSTER.cpp -o debug_sigmerIndex_to_binary $(LDFLAGS)


clean:
	rm -rf *.o  sigmer_generation sigmer_count_01_new EM
