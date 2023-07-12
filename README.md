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
