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
#include <gmpxx.h>
#include <chrono>
#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

using namespace seqan2;
using namespace std;
using namespace std::chrono;

typedef std::pair<CharString, int> location;
typedef std::size_t position_t;

namespace bip = boost::interprocess;

typedef bip::allocator<unsigned long long, bip::managed_mapped_file::segment_manager> ShmemAllocator;
typedef bip::vector<unsigned long long, ShmemAllocator> MyVector;

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

int encodeGMP(CharString &kmer, mpz_class &x_mpz)
{
    for (size_t i = 0; i < length(kmer); ++i) {
	if(kmer[i] == 'A')
	{
		mpz_clrbit(x_mpz.get_mpz_t(), 2*i); // false
		mpz_clrbit(x_mpz.get_mpz_t(), (2*i)+1); // false		
	}
	else if(kmer[i] == 'G')
	{
		mpz_clrbit(x_mpz.get_mpz_t(), 2*i); // false
		mpz_setbit(x_mpz.get_mpz_t(), (2*i)+1); // true
	}
	else if(kmer[i] == 'C')
	{
		mpz_setbit(x_mpz.get_mpz_t(), 2*i); // true
		mpz_clrbit(x_mpz.get_mpz_t(), (2*i)+1); // false			
	}
	else if(kmer[i] == 'T')
	{
		mpz_setbit(x_mpz.get_mpz_t(), 2*i); // true
		mpz_setbit(x_mpz.get_mpz_t(), (2*i)+1); // true			    
	}
	else
	{
		return 1;
	}
    }
    return 0;
}

int encode(const std::string& kmer, unsigned long long int& x_ull) 
{
	x_ull = 0;
	for (size_t i = 0; i < kmer.length(); ++i) 
	{
		char nucleotide = kmer[i];
		switch (nucleotide) 
		{
			case 'A':
				// No need to change bits, as both bits are 0 by default in an unsigned long long int
				break;
			case 'G':
				x_ull |= (1ULL << (2*i + 1));
				break;
			case 'C':
				x_ull |= (1ULL << (2*i));
				break;
			case 'T':
				x_ull |= (1ULL << (2*i));
				x_ull |= (1ULL << (2*i + 1));
				break;
			default:
				return 1; // Return an error code if an invalid character is encountered
		}
	}
	return 0;
}

int readFile(const char* filename)
{
	//Read numbers from the binary file and check if they match the originals
	FILE* file = fopen(filename, "rb"); // Open file in binary read mode
	if (file == nullptr) 
	{
		std::cerr << "Error opening file for reading." << std::endl;
		return 1;
	}

	mpz_t number;
	mpz_init(number);

	while (!feof(file))
	{
		size_t bytes_read = mpz_inp_raw(number, file);
		if (bytes_read == 0) 
		{
			if (feof(file))
			{
				// end of file reached
				break;
			}
			else
			{
				cerr << "ERROR reading from file." << endl;
				fclose(file);
				mpz_clear(number);
				return 1;
			}
		}	
		//gmp_printf("Read number: %Zd\n", number);
	}
	fclose(file);
	mpz_clear(number);
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
		cout << "Reading kmer db " << i << endl;
  		KMC::Runner runner;
  		CKMCFile kmer_database;

		auto start = high_resolution_clock::now();

		bip::managed_mapped_file file(bip::open_or_create, "shared_memory.bin", 10); // Adjust size as needed
		// Construct the allocator
		ShmemAllocator alloc(file.get_segment_manager());
		// Construct the vector in shared memory
		MyVector *vec = file.find_or_construct<MyVector>("Vector")(alloc);

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

		// each kmer
  		while(kmer_database.ReadNextKmer(kmer_object, counter))
  		{
			string str;

			if(counter >= 5 && counter <= 15)
			{	
				kmer_object.to_string(str);
				unsigned long long int encoded;
				encode(str, encoded);
				//kmers.push_back(encoded);
				vec->push_back(encoded);
				cout << "encoded\t" << encoded << endl; 
			}
	  	}

		for(const auto& n : *vec)
		{
			cout << "unsorted\t" << n << endl;
		}

		kmer_database.Close();
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
		cout << "File read: " << duration.count() << endl;

		// sort our vector
		start = high_resolution_clock::now();
		std::sort(vec->begin(), vec->end());
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		cout << "Sorted : " << duration.count() << endl;

		for(const auto& n : *vec)
		{
			cout << "done_sorted\t" << n << endl;
		}
	}	

	return 0;

}
