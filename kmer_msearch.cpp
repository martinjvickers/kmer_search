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


using namespace seqan2;
using namespace std;

typedef std::pair<CharString, int> location;
typedef std::size_t position_t;

struct ModifyStringOptions {
  CharString matrixFileName;
  CharString kmer;
  CharString seqFile;
  vector<CharString> databases;
  int threads;
};

seqan2::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                                    int argc, char const ** argv) {
  ArgumentParser parser("kmer_msearch");

  addOption(parser, ArgParseOption("s", "seq-file",
                                   "Path to the input file containing sequences",
                                   ArgParseArgument::INPUT_FILE,
                                   "IN"));
  setRequired(parser, "seq-file");
/*
  addOption(parser, ArgParseOption("db", "databases",
  	  			   "A list of databases.",
				   ArgParseArgument::STRING, "TEXT", true));
  setRequired(parser, "databases");
*/

  addOption(parser, ArgParseOption("t", "threads", "Number of threads",
                                   ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "threads", "1");

  setShortDescription(parser, "kmer_msearch");
  setVersion(parser, "0.0.1");
  setDate(parser, "June 2023");
  addUsageLine(parser, "-db kmerdb1 -db kmerdb2 -s seq.fa \
                        [\\fIOPTIONS\\fP] ");
  addDescription(parser, "");
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  if (res != ArgumentParser::PARSE_OK){
    return res;
  }

  getOptionValue(options.seqFile, parser, "seq-file");
/*  for(int i = 0; i < getOptionValueCount(parser, "databases"); i++)
  {
    CharString tmpVal;
    getOptionValue(tmpVal, parser, "databases", i);
    options.databases.push_back(tmpVal);
  }	   
*/
  return ArgumentParser::PARSE_OK;
}

class CountMinSketch {
private:
    std::vector<std::vector<uint64_t>> sketch;
    std::vector<uint64_t> hash_coefficients;
    size_t d; // Depth (number of rows)
    size_t w; // Width (number of counters per row)

public:
    CountMinSketch(size_t depth, size_t width)
        : d(depth), w(width), sketch(depth, std::vector<uint64_t>(width, 0)) {
        // Initialize the hash coefficients
        std::hash<uint64_t> hash_func;
        for (size_t i = 0; i < d; ++i) {
            uint64_t coefficient = hash_func(i); // Use different hash coefficients for each row
            hash_coefficients.push_back(coefficient);
        }
    }

    // Increment the count for the given value
    void increment(uint64_t value) {
        for (size_t i = 0; i < d; ++i) {
            // Compute the hash index for the current row
            uint64_t hash_index = (hash_coefficients[i] ^ value) % w;
            sketch[i][hash_index]++;
        }
    }

    // Estimate the count for the given value
    uint64_t estimate(uint64_t value) const {
        uint64_t min_count = UINT64_MAX;
        for (size_t i = 0; i < d; ++i) {
            // Compute the hash index for the current row
            uint64_t hash_index = (hash_coefficients[i] ^ value) % w;
            min_count = std::min(min_count, sketch[i][hash_index]);
        }
        return min_count;
    }
};

/*
class CountMinSketch {
private:
    std::vector<std::vector<uint64_t>> sketch;
    std::vector<uint64_t> hash_coefficients;
    size_t d; // Depth (number of rows)
    size_t w; // Width (number of counters per row)
    size_t num_threads; // Number of threads to use for parallelization

    // Helper struct to pass data to threads
    struct ThreadData {
        CountMinSketch* cms;
        uint64_t value;
        size_t thread_id;
    };

    // Helper function to increment a subset of rows
    static void* incrementSubset(void* arg) {
        ThreadData* data = static_cast<ThreadData*>(arg);
        size_t rows_per_thread = data->cms->d / data->cms->num_threads;
        size_t remaining_rows = data->cms->d % data->cms->num_threads;
        size_t start_row = data->thread_id * rows_per_thread;
        size_t end_row = start_row + rows_per_thread + (data->thread_id == data->cms->num_threads - 1 ? remaining_rows : 0);

        for (size_t i = start_row; i < end_row; ++i) {
            uint64_t hash_index = (data->cms->hash_coefficients[i] ^ data->value) % data->cms->w;
            data->cms->sketch[i][hash_index]++;
	    cout << i << "\t" << hash_index << "\t" << data->cms->sketch[i][hash_index] << endl;
        }

        return nullptr;
    }

public:
    CountMinSketch(size_t depth, size_t width, size_t num_threads)
        : d(depth), w(width), num_threads(num_threads),
          sketch(depth, std::vector<uint64_t>(width, 0)) {
        // Initialize the hash coefficients
        std::hash<uint64_t> hash_func;
        for (size_t i = 0; i < d; ++i) {
            uint64_t coefficient = hash_func(i); // Use different hash coefficients for each row
            hash_coefficients.push_back(coefficient);
        }
    }

    // Increment the count for the given value (parallelized version)
    void increment(uint64_t value) {
        std::vector<pthread_t> threads(num_threads);
        std::vector<ThreadData> thread_data(num_threads);

        cout << num_threads << endl;

        // Launch threads
        for (size_t i = 0; i < num_threads; ++i) {
            thread_data[i] = {this, value, i};
	    cout << value << "\t" << i << "\t" << num_threads << endl;
            pthread_create(&threads[i], nullptr, incrementSubset, &thread_data[i]);
        }

        // Wait for threads to finish
        for (size_t i = 0; i < num_threads; ++i) {
            pthread_join(threads[i], nullptr);
        }
    }

    // Estimate the count for the given value
    uint64_t estimate(uint64_t value) const {
      uint64_t min_count = UINT64_MAX;
      for (size_t i = 0; i < d; ++i) {
        // Compute the hash index for the current row
        uint64_t hash_index = (hash_coefficients[i] ^ value) % w;
        min_count = std::min(min_count, sketch[i][hash_index]);
      }
      return min_count;
    }

    void printSketch(){
      for(size_t i = 0; i < d; ++i) 
      {
        for(size_t j = 0; j < w; ++j)
	{
          cout << sketch[i][j] << "\t";		
        }
	cout << endl;
      }	       
    }
};
*/

int encode(CharString &kmer, bitset<62> &kmer2)
{
  for(int i=0; i < length(kmer); i++)
  { 
    if(kmer[i] == 'A')
    {
      kmer2.set(position_t(2*i), false);
      kmer2.set(position_t((2*i)+1), false);
    }
    else if(kmer[i] == 'G')
    {
      kmer2.set(position_t(2*i), false);
      kmer2.set(position_t((2*i)+1), true);
    }
    else if(kmer[i] == 'C')
    {
      kmer2.set(position_t(2*i), true);
      kmer2.set(position_t((2*i)+1), false);
    }
    else if(kmer[i] == 'T')
    {
      kmer2.set(position_t(2*i), true);
      kmer2.set(position_t((2*i)+1), true);
    }
    else
    {
      return 1;
    }
  }
  return 0;
}

unsigned long long int ipow( unsigned long long int base, unsigned long long int exp )
{
  unsigned long long int result = 1ULL;
  while( exp )
  {
    if ( exp & 1 )
    {
      result *= (unsigned long long int)base;
    }
    exp >>= 1;
    base *= base;
  }
    return result;
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

  int _kmer_length = 31;
  SeqFileIn seqFileIn(toCString(options.seqFile));
  int read_count = 0;
  long gb = 1024 * 1024 * 1024;

  //CountMinSketch cms(1000, 1024, options.threads);
  CountMinSketch cms(1000, 1024);

  while (!atEnd(seqFileIn))
  {
    CharString id, inf;
    Dna5String seq;
    readRecord(id, seq, seqFileIn);

    //cout << seq << endl;

    for(int i = 0; i < length(seq)-_kmer_length+1; i++)
    {
      inf = infix(seq, i, i+_kmer_length);
      std::bitset<62> kmer;
      encode(inf, kmer);
      //cout << kmer.to_ullong() << endl;           

      cms.increment((uint64_t)kmer.to_ullong());
      //cout << inf << "\t" << kmer.to_ulong() << endl;
      //cout << hash_value(kmer) << endl;
    }
  }

  //cms.printSketch();

  close(seqFileIn);

  SeqFileIn seqFileIn2(toCString(options.seqFile));

  while (!atEnd(seqFileIn2))
  {
    CharString id, inf;
    Dna5String seq;
    readRecord(id, seq, seqFileIn2);

    for(int i = 0; i < length(seq)-_kmer_length+1; i++)
    {
      inf = infix(seq, i, i+_kmer_length);
      std::bitset<62> kmer;
      encode(inf, kmer);
      cout << inf << "\t" << kmer.to_ullong() << "\t" << cms.estimate((uint64_t)kmer.to_ullong()) << endl;
    }    
  }
  close(seqFileIn2);
  return 0;

}
