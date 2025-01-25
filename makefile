#CC = g++ -O3 -march=native -std=c++14 -mavx2 -DSEQAN_HAS_ZLIB -DNDEBUG -DSEQAN_DISABLE_VERSION_CHECK -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -static -static-libgcc -static-libstdc++
#CC = g++ -O3 -march=native -std=c++14 -DSEQAN_HAS_ZLIB -DNDEBUG -DSEQAN_DISABLE_VERSION_CHECK -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -static -static-libgcc -static-libstdc++
CC = g++ -O3 -std=c++14 -DSEQAN_HAS_ZLIB -DNDEBUG -DSEQAN_DISABLE_VERSION_CHECK -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -static -static-libgcc -static-libstdc++
LIB = -Wl,--whole-archive -lpthread -lgmp -lgmpxx -lz -lm -llzma -lbz2 -lrt -Wl,--no-whole-archive
INC = -I./seqan/include
LD =  KMC/bin/libkmc_core.a 

all : KMC/bin/libkmc_core.a kmer_search kmer_pa_matrix kmer_pa_matrix_new kmer_pa_matrix_search kmer_pa_matrix_gene_search kmer_pa_matrix_dump kmer_pa_matrix_gene_search_multi

kmer_search : kmer_search.cpp
		$(CC) -O3 $^ $(LD) $(INC) $(LIB) -o $@
kmer_pa_matrix : kmer_pa_matrix.cpp
		$(CC) -O3 $^ $(LD) $(INC) $(LIB) -o $@
#kmer_msearch : kmer_msearch.cpp MurmurHash3.h  MurmurHash3.cpp
#		$(CC) -O3 $^ $(LD) $(INC) $(LIB) -o $@
kmer_pa_matrix_new : kmer_pa_matrix_new.cpp
		$(CC) -O3 $^ $(LD) $(INC) $(LIB) -o $@
kmer_pa_matrix_search : kmer_pa_matrix_search.cpp
		$(CC) -O3 $^ $(LD) $(INC) $(LIB) -o $@
kmer_pa_matrix_gene_search : kmer_pa_matrix_gene_search.cpp
		$(CC) -O3 $^ $(LD) $(INC) $(LIB) -o $@
kmer_pa_matrix_dump : kmer_pa_matrix_dump.cpp
		$(CC) -O3 $^ $(LD) $(INC) $(LIB) -o $@
kmer_pa_matrix_gene_search_multi : kmer_pa_matrix_gene_search_multi.cpp
		$(CC) -O3 $^ $(LD) $(INC) $(LIB) -o $@
clean:
		rm kmer_search kmer_msearch kmer_pa_matrix kmer_pa_matrix_new kmer_pa_matrix_search kmer_pa_matrix_gene_search kmer_pa_matrix_dump kmer_pa_matrix_gene_search_multi
