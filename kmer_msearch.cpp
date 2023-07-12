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


using namespace seqan2;
using namespace std;

typedef std::pair<CharString, int> location;
typedef std::size_t position_t;

struct ModifyStringOptions {
  CharString matrixFileName;
  CharString kmer;
  CharString seqFile;
  vector<CharString> databases;
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

int encode(CharString &kmer, bitset<16> &kmer2)
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

  // let's have some 8-mer arrays, four of them to deal with 31-mers
  int* arrayA = (int*)malloc(pow(2, 8) * sizeof(int));
  memset(arrayA, 0, pow(2, 8) * sizeof(int));
  int* arrayB = (int*)malloc(pow(2, 8) * sizeof(int));
  memset(arrayB, 0, pow(2, 8) * sizeof(int));
  int* arrayC = (int*)malloc(pow(2, 8) * sizeof(int));
  memset(arrayC, 0, pow(2, 8) * sizeof(int));
  int* arrayD = (int*)malloc(pow(2, 8) * sizeof(int));
  memset(arrayD, 0, pow(2, 8) * sizeof(int));

  //cout << (pow(2, (_kmer_length*2)) * sizeof(int)) / gb << " ram in GB" << endl;
  //cout << (ipow(2, (_kmer_length*2)) * sizeof(int)) / gb << " ram in GB" << endl;

  //cout << "allocating " << pow(2, _kmer_length) << endl;
  //int* array = (int*)malloc(pow(2, _kmer_length) * sizeof(int));
  //memset(array, 0, pow(2, _kmer_length) * sizeof(int));
  //cout << "array " << sizeof(array)/sizeof(array[0]) << endl; 
 
/*  
  unsigned long long int combinations = ipow(2, (_kmer_length*2));
  unsigned short int* array;
  cout << "number of elements " << combinations << endl;
  cout << combinations * sizeof(unsigned short int) /gb << " GB" << endl;
  cout << "bool size in bytes " << sizeof(unsigned short int) << endl;

  try {
    array = new unsigned short int[combinations];
    memset(array, 0, combinations * sizeof(unsigned short int));
    cout << "Successfully allocated memory" << endl;
    //cout << "combinations=" << combinations << "\tcombinations-1=" << array[combinations-1] << endl;
  }
  catch(const bad_alloc& e) {
    cout << "alloc failed : " << e.what() << endl;
    return 1;    
  }	  
  cout << "allocated" << endl;
  */

/*
  for(unsigned long long int i = 0; i < combinations; i++)
  {
    cout << i << "\t" << array[i] << endl;	  
  }	  
*/

  while (!atEnd(seqFileIn))
  {
    CharString id;
    Dna5String seq;
    readRecord(id, seq, seqFileIn);
    //cout << seq << endl;
    CharString inf;

    for(int i = 0; i < length(seq)-_kmer_length+1; i++)
    {
      inf = infix(seq, i, i+7);
      std::bitset<16> kmer2;

      if(encode(inf, kmer2) == 0)
      {
        cout << kmer2.to_ulong() << endl;	      
	arrayA[kmer2.to_ulong()]++;	
      }

      inf = infix(seq, i+8, i+16);
      if(encode(inf, kmer2) == 0)
      { 
	 cout << kmer2.to_ulong() << endl;     
         arrayB[kmer2.to_ulong()]++;
      }
     
    }
    read_count++;
  }

  cout << "completed " << read_count << " reads." << endl;

/*  
  for(auto d : options.databases)
  {	  

    KMC::Runner runner;
    CKMCFile kmer_database;

    if(!kmer_database.OpenForRA(toCString(d)))
    {
      cout <<"Error opening kmer database file." << endl;
    }

    uint32 _kmer_length, _mode, _counter_size, _lut_prefix_length, _sig_len, _min_count;
    uint64 _max_count;
    uint64 _total_kmers;

    kmer_database.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _sig_len, _min_count, _max_count, _total_kmers);

    SeqFileIn seqFileIn(toCString(options.seqFile));

    while (!atEnd(seqFileIn))
    {
      CharString id;
      Dna5String seq;
      readRecord(id, seq, seqFileIn);
      CharString inf;

      for(int i = 0; i < length(seq)-_kmer_length+1; i++)
      {
         inf = infix(seq, i, i+_kmer_length);

         uint64 counter;

         CKmerAPI kmer_object(_kmer_length);
         std::string tofind = toCString(inf);
         kmer_object.from_string(tofind);

         if(kmer_database.CheckKmer(kmer_object, counter))
         {
           cout << kmer_object.to_string() << "\t" << counter << endl;
         }
         else
         {
           kmer_object.reverse();
           if(kmer_database.CheckKmer(kmer_object, counter))
           {
             cout << kmer_object.to_string() << "\t" << counter << endl;
           }
         }
      }
    }
    kmer_database.Close();
  }
*/
//  delete [] array;
//
  delete [] arrayA;
  delete [] arrayB;
  delete [] arrayC;
  delete [] arrayD;
  return 0;

}
