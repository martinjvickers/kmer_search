/*
 *  This is not part of the program, this is me faffing around 
 *  with things
 *
 *  Martin Vickers
 *
 *
 * */

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
#include <seqan/gff_io.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <boost/algorithm/string.hpp>
#include <map>
#include <utility>
#include <iterator>
#include <random>
#include <algorithm>
#include "KMC/include/kmc_runner.h"
#include "KMC/kmc_api/kmc_file.h"
#include <cmath>
#include <bitset>
#include <functional>
#include <pthread.h>
#include "MurmurHash3.h"
//#include <immintrin.h>
#include <chrono>

using namespace std::chrono;
using namespace seqan2;
using namespace std;

typedef std::pair<CharString, int> location;
typedef std::size_t position_t;

struct ModifyStringOptions {
  CharString matrixFileName;
  CharString kmer;
  CharString seqFile;
  vector<CharString> databases;
  int depth;
  int width;
};

seqan2::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                                    int argc, char const ** argv) {
  ArgumentParser parser("kmer_msearch");

  addOption(parser, ArgParseOption("s", "seq-file",
                                   "Path to the input file containing sequences",
                                   ArgParseArgument::INPUT_FILE,
                                   "IN"));
  setRequired(parser, "seq-file");
  addOption(parser, ArgParseOption("w", "width", "width",
                                   ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "width", "1024");
  addOption(parser, ArgParseOption("d", "depth", "depth",
                                   ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "depth", "10000");

  setShortDescription(parser, "kmer_msearch");
  setVersion(parser, "0.0.1");
  setDate(parser, "October 2023");
  addUsageLine(parser, "-s seq.fa -w 16384 -d 1000000 \
                        [\\fIOPTIONS\\fP] ");
  addDescription(parser, "");
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  if (res != ArgumentParser::PARSE_OK){
    return res;
  }

  getOptionValue(options.seqFile, parser, "seq-file");
  getOptionValue(options.depth, parser, "depth");
  getOptionValue(options.width, parser, "width");
  return ArgumentParser::PARSE_OK;
}

void increment(ModifyStringOptions options, unsigned long long value, uint16_t** sketch,
               std::vector<uint16_t> &hash_coefficients)
{
  for (int i = 0; i < options.depth; ++i)
  {
    // Compute the hash index for the current row
    uint16_t hash_index = (hash_coefficients[i] ^ value) % options.width;
    sketch[i][hash_index]++;
  }
}

// Estimate the count for the given value
uint16_t estimate(ModifyStringOptions options, unsigned long long value, uint16_t** sketch,
		  std::vector<uint16_t> &hash_coefficients)
{
  uint16_t min_count = UINT16_MAX;
  for (int i = 0; i < options.depth; ++i)
  {
    // Compute the hash index for the current row
    uint16_t hash_index = (hash_coefficients[i] ^ value) % options.width;
    min_count = std::min(min_count, sketch[i][hash_index]);
  }
  return min_count;
}

int encode(CharString &kmer, bitset<102> &kmer2)
{
  for(int i=0; i < length(kmer); i++)
  { 
    switch(kmer[i])
    {
      case 'A':
        kmer2.set(position_t(2*i), false);
	kmer2.set(position_t((2*i)+1), false);
      case 'G':
	kmer2.set(position_t(2*i), false);
	kmer2.set(position_t((2*i)+1), true);
      case 'C':
        kmer2.set(position_t(2*i), true);
        kmer2.set(position_t((2*i)+1), false);
      case 'T':
        kmer2.set(position_t(2*i), true);
        kmer2.set(position_t((2*i)+1), true);
      default:
        return 1;
    }    
  }
  return 0;
}

// A basic template to get up and running quickly
int main(int argc, char const ** argv)
{
  //parse our options
  ModifyStringOptions options;
  ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
  if(res != ArgumentParser::PARSE_OK){
    return res == ArgumentParser::PARSE_ERROR;
  }

  int _kmer_length = 51;
  SeqFileIn seqFileIn(toCString(options.seqFile));
  int read_count = 0;

  map<uint64_t, uint16_t> store;
  auto start = high_resolution_clock::now(); // timer

  //set up the sketch
  uint16_t** sketch = new uint16_t*[options.depth];
  for (size_t i = 0; i < options.depth; ++i) {
    sketch[i] = new uint16_t[options.width];
    memset(sketch[i], 0, options.width * sizeof(uint16_t));
  }

  //set up the coefficients
  std::vector<uint16_t> hash_coefficients;
  for (int i = 0; i < options.depth; ++i) {
    uint64_t coefficient;
    MurmurHash3_x64_128(&i, sizeof(i), i, &coefficient);
    hash_coefficients.push_back(coefficient);
  }

  // count each kmer
  while (!atEnd(seqFileIn))
  {
    CharString id, inf;
    Dna5String seq;
    readRecord(id, seq, seqFileIn);

    for(int i = 0; i < length(seq)-_kmer_length+1; i++)
    {
      inf = infix(seq, i, i+_kmer_length);
      std::bitset<102> kmer;
      encode(inf, kmer);

      //auto start = high_resolution_clock::now();
      increment(options, kmer.to_ullong(), sketch, hash_coefficients);
      //auto stop = high_resolution_clock::now();
      //auto duration = duration_cast<microseconds>(stop - start);
      store[kmer.to_ullong()]++;
    }
  }

  close(seqFileIn);

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "Finished counting - time=" << duration.count() << "ms" << endl;
  cout << "Calculating accuracy" << endl;

  // this is if we want to calculate our accuracy
  SeqFileIn seqFileIn2(toCString(options.seqFile));
  size_t num_correct_estimations = 0;
  size_t num_elements = 0;

  while (!atEnd(seqFileIn2))
  {
    CharString id, inf;
    Dna5String seq;
    readRecord(id, seq, seqFileIn2);

    for(int i = 0; i < length(seq)-_kmer_length+1; i++)
    {
      inf = infix(seq, i, i+_kmer_length);
      std::bitset<102> kmer;
      encode(inf, kmer);
      //cout << inf << "\t" << kmer.to_ullong() << "\t" << estimate(options, kmer.to_ullong(), sketch, hash_coefficients) << endl;
      //
      if(estimate(options, kmer.to_ullong(), sketch, hash_coefficients) == store[kmer.to_ullong()])
      {
        num_correct_estimations++;
      }
      num_elements++;
    }    
  }

  double accuracy = static_cast<double>(num_correct_estimations) / num_elements * 100.0;
  cout << "Number of k-mers\t" << num_elements << endl;
  cout << "Number of correct counts\t" << num_correct_estimations << endl;
  cout << "Accuracy\t" << accuracy << endl;

  close(seqFileIn2);
  // end accuracy


  // clean up
  for (size_t i = 0; i < options.depth; ++i) {
    delete[] sketch[i];
  }
  delete[] sketch;
  return 0;

}
