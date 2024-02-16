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
#include <boost/dynamic_bitset.hpp>

using namespace seqan2;
using namespace std;

typedef std::pair<CharString, int> location;
typedef std::size_t position_t;

struct ModifyStringOptions {
  vector<CharString> kmerDatabases;
};

seqan2::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                                    int argc, char const ** argv) {
  ArgumentParser parser("kmer_pa_matrix");
  addOption(parser, ArgParseOption("k", "kmer-database",
			           "A list of kmer databases",
                                    ArgParseArgument::STRING, "TEXT", true));
  setRequired(parser, "kmer-database");


  setShortDescription(parser, "kmer_pa_matrix");
  setVersion(parser, "0.0.1");
  setDate(parser, "January 2024");
  addUsageLine(parser, "-k kmerdb \
                        [\\fIOPTIONS\\fP] ");
  addDescription(parser, "");
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  if (res != ArgumentParser::PARSE_OK){
    return res;
  }

  for(int i = 0; i < getOptionValueCount(parser, "kmer-database"); i++)
  {
	CharString tmpVal;
	getOptionValue(tmpVal, parser, "kmer-database", i);
	options.kmerDatabases.push_back(tmpVal);
  }

  return ArgumentParser::PARSE_OK;
}

int encode(CharString &kmer, boost::dynamic_bitset<> &kmer2)
{
	for(int i = 0; i < length(kmer); i++)
	{
	/*
	 * I don't know why switches aren't working, i will deal with this later
	 * switch(kmer[i])
		{
		  case 'A':
		    cout << "A\t";
		    kmer2[2*i] = false;
                    kmer2[(2*i)+1] = false;
		  case 'G':
		    cout << "G\t";
                    kmer2[2*i] = false;
                    kmer2[(2*i)+1] = true;
                  case 'C':
		    cout << "C\t";
	            kmer2[2*i] = true;
                    kmer2[(2*i)+1] = false;
                  case 'T':
		    cout << "T\t";
                    kmer2[2*i] = true;
                    kmer2[(2*i)+1] = true;
                }
        }
	*/
		if(kmer[i] == 'A')
		{
			kmer2[2*i] = false;
                        kmer2[(2*i)+1] = false;
		}
		else if(kmer[i] == 'G')
		{
                    kmer2[2*i] = false;
                    kmer2[(2*i)+1] = true;		    
		}
		else if(kmer[i] == 'C')
		{
                    kmer2[2*i] = true;
                    kmer2[(2*i)+1] = false;		    
		}
		else if(kmer[i] == 'T')
		{
                    kmer2[2*i] = true;
                    kmer2[(2*i)+1] = true;		    
		}
		else
		{
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

  // ideally, this loop would be scaled across multiple threads/machines
  for(auto i : options.kmerDatabases)
  {	  
  	KMC::Runner runner;
  	CKMCFile kmer_database;

  	if(!kmer_database.OpenForListing(toCString(i)))
  	{
  		cout <<"Error opening kmer database file." << endl;
  	} 

	// read existing kmer database, which gives us our kmer length
  	uint32 _kmer_length, _mode, _counter_size, _lut_prefix_length, _sig_len, _min_count;
  	uint64 _max_count;
  	uint64 _total_kmers;

  	kmer_database.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _sig_len, _min_count, _max_count, _total_kmers);

  	CKmerAPI kmer_object(_kmer_length);
  	uint32 counter;

        int limit = 10;
        int c = 0;

	std::map<std::vector<bool>, bool> store_kmers;

	// each kmer
  	while(kmer_database.ReadNextKmer(kmer_object, counter))
  	{
		string str;

		if(counter > 2 && counter < 15)
		{	
			kmer_object.to_string(str);

			CharString kmer = str;
			boost::dynamic_bitset<> x(_kmer_length*2);
			if(encode(kmer, x) == 1)
			{
				cout << "Encoding of " << kmer << " failed" << endl;
			}
		}
  	}		  
	kmer_database.Close();
  }	
  return 0;

}
