/*
 *  kmer_search
 *  Written by Martin Vickers @ John Innes Centre
 *  martin.vickers @ jic dot ac dot uk
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

using namespace seqan2;
using namespace std;

typedef std::pair<CharString, int> location;

struct ModifyStringOptions {
  CharString kmerDatabase;
  CharString seqFile;
};

seqan2::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                                    int argc, char const ** argv) {
  ArgumentParser parser("kmer_search");
  addOption(parser, ArgParseOption("k", "kmer-database", 
                                   "Path to the input file", 
                                   ArgParseArgument::INPUT_FILE,
                                   "IN"));
  setRequired(parser, "kmer-database");

  addOption(parser, ArgParseOption("f", "fasta-file",
                                   "Path to the input file containing sequences",
                                   ArgParseArgument::INPUT_FILE,
                                   "IN"));
  setRequired(parser, "fasta-file");

  setShortDescription(parser, "kmer_search");
  setVersion(parser, "0.0.1");
  setDate(parser, "July 2023");
  addUsageLine(parser, "-k kmerdb -f seq.fa \
                        [\\fIOPTIONS\\fP] ");
  addDescription(parser, "");
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  if (res != ArgumentParser::PARSE_OK){
    return res;
  }

  getOptionValue(options.kmerDatabase, parser, "kmer-database");
  getOptionValue(options.seqFile, parser, "fasta-file");

  return ArgumentParser::PARSE_OK;
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

  KMC::Runner runner;
  CKMCFile kmer_database;

  if(!kmer_database.OpenForRA(toCString(options.kmerDatabase)))
  {
    cout <<"Error opening kmer database file." << endl;
  } 

  // read existing kmer database, which gives us our kmer length
  uint32 _kmer_length, _mode, _counter_size, _lut_prefix_length, _sig_len, _min_count;
  uint64 _max_count;
  uint64 _total_kmers;

  kmer_database.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _sig_len, _min_count, _max_count, _total_kmers);

  SeqFileIn seqFileIn(toCString(options.seqFile));

  /* loop through fasta file, extract kmers of kmer_length, print out the occurances of
   * said kmers in the database(s)
   *
   */
  while (!atEnd(seqFileIn))
  {
    CharString id;
    Dna5String seq;
    readRecord(id, seq, seqFileIn);
    CharString inf;

    for(int i = 0; i < length(seq)-_kmer_length+1; i++)
    {
       inf = infix(seq, i, i+_kmer_length);
       //cout << inf << endl;
       //cout << inf << "\t" << length(inf) << "\t" << i << "\t" << i+_kmer_length << "\t" << _kmer_length << endl;

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
         //kmer not in the database
       }
    }
  }

  return 0;

}
