# kmer_search

### Compiling

Getting this to compile on my machine;

```
git submodule update --init --recursive
cd KMC
make -j8
cd ../
make
```

### Running

```
./kmer_search -f examples/test_dna.fq -k examples/3mer
```




```
singularity exec ~/oldscratch/bin/kmer_pa.img kmer_pa_matrix_new -k meh.txt -m master_kmers > result.txt
singularity exec ~/oldscratch/bin/kmer_pa.img kmer_pa_matrix_new -k meh.txt -m master_kmers > result.txt
```



```
kmer_pa_matrix_new -k examples/massive_list.txt -m examples/3mers -o meh.out --create-pa
kmer_pa_matrix_new -k examples/massive_list.txt -m examples/3mers -o meh.out -n NUMBER_OF_RECORD
```
