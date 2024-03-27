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
#include "KMC/include/kmc_runner.h"
#include "KMC/kmc_api/kmc_file.h"
#include <chrono>

using namespace seqan2;
using namespace std;
using namespace std::chrono;

struct ModifyStringOptions 
{
	CharString kmerDatabasesFilenames;
	CharString masterKmerDatabase;
};

struct kmer
{
	uint64 k;
	unsigned short bits;
};

seqan2::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                                    int argc, char const ** argv) {
  ArgumentParser parser("kmer_pa_matrix");
  addOption(parser, ArgParseOption("k", "kmer-databases",
			           "A file with a list of kmer databases",
                                    ArgParseArgument::STRING, "TEXT", true));
  setRequired(parser, "kmer-databases");

  addOption(parser, ArgParseOption("m", "master-database",
			           "A master list of the union of all kmer databases",
                                   ArgParseArgument::STRING, "TEXT", true));
  setRequired(parser, "master-database");

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

  getOptionValue(options.kmerDatabasesFilenames, parser, "kmer-databases");
  getOptionValue(options.masterKmerDatabase, parser, "master-database");

  return ArgumentParser::PARSE_OK;
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

int decode(const unsigned long long int& x_ull, std::string& kmer, uint32 kmer_length) 
{
	kmer.clear();
	for (size_t i = 0; i < kmer_length; ++i)
	{
		unsigned long long int mask = 3ULL << (2*i);
		unsigned long long int nucleotide_bits = (x_ull & mask) >> (2*i);
		switch (nucleotide_bits)
		{
			case 0:
				kmer += 'A';
				break;
			case 1:
				kmer += 'C';
				break;
			case 2:
				kmer += 'G';
				break;
			case 3:
				kmer += 'T';
				break;
			default:
				return 1; // Return an error code if unexpected bits are encountered
		}
	}
	return 0;
}

int checkForPA(CharString masterKmerDB, CharString currentDBtoCheck, vector<kmer> &v, int n)
{
	// open master kmer file in dump mode
	// open kmer db we are checking in random mode
	// loop through master in order, checking to see if it exists
	// flip the bit
	KMC::Runner runner;
	CKMCFile kmer_master_database;
	if(!kmer_master_database.OpenForListing(toCString(masterKmerDB)))
	{
		cout <<"Error opening kmer database file " << masterKmerDB<< endl;
		return 1;
	}
	uint32 _kmer_length, _mode, _counter_size, _lut_prefix_length, _sig_len, _min_count;
	uint64 _max_count, _total_kmers, counter;
	kmer_master_database.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _sig_len, _min_count, _max_count, _total_kmers);
	CKmerAPI kmer_object(_kmer_length);

	// this is the query database for random access searching
	CKMCFile kmer_query_database;
	if(!kmer_query_database.OpenForRA(toCString(currentDBtoCheck)))
	{
		cout <<"Error opening kmer database file " << currentDBtoCheck << endl;
		return 1;
	}
	uint32 _kmer_query_length, _query_mode, _query_counter_size, _query_lut_prefix_length, _query_sig_len, _query_min_count;
	uint64 _query_max_count, _query_total_kmers;
	kmer_query_database.Info(_kmer_query_length, _query_mode, _query_counter_size, _query_lut_prefix_length, _query_sig_len, _query_min_count, _query_max_count, _query_total_kmers);
	CKmerAPI query_kmer_object(_kmer_query_length);

	int c = 0;

	// loop through each kmer
	while(kmer_master_database.ReadNextKmer(kmer_object, counter))
	{
		if(kmer_query_database.IsKmer(kmer_object) == 1)
		{
			v[c].bits ^= (static_cast<uint64_t>(1) << n);
		}
		c++;
	}

	kmer_master_database.Close();
	kmer_query_database.Close();
	return 0;
}

vector<CharString> createFileList(CharString kmerDatabaseList)
{
	vector<CharString> v;
	string line;
	ifstream file(toCString(kmerDatabaseList));
	if (file.is_open())
	{
		while (getline(file, line)) 
		{
			v.push_back(line);
		}
		file.close();
	}
	return v;
}

int createOutputArray(CharString masterKmers, int n, vector<kmer> &v, uint32 &kmer_length)
{
	auto start = high_resolution_clock::now();
	KMC::Runner runner;
	CKMCFile kmer_database;
	if(!kmer_database.OpenForListing(toCString(masterKmers)))
	{
		cout <<"Error opening kmer database file." << endl;
		return 1;
	}
	// read existing kmer database, which gives us our kmer length
	uint32 _kmer_length, _mode, _counter_size, _lut_prefix_length, _sig_len, _min_count;
	uint64 _max_count;
	uint64 _total_kmers;
	kmer_database.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _sig_len, _min_count, _max_count, _total_kmers);
	kmer_length = _kmer_length;
	CKmerAPI kmer_object(_kmer_length);
	uint32 counter;

	cerr << "Master database consists of " << _total_kmers << " " << _kmer_length << "-mers"<< endl;
	v.resize(_total_kmers);
	uint64 c = 0;

	while(kmer_database.ReadNextKmer(kmer_object, counter))
	{
		string str;
		kmer_object.to_string(str);
		unsigned long long int encoded;
		encode(str, encoded);
		//kmer kmerToAdd;
		//kmerToAdd.k = encoded;
		//unsigned short zeroed = 0;
		//kmerToAdd.bits = zeroed;
		//v.push_back(kmerToAdd);
		v[c].k = encoded;
		v[c].bits = 0;
		c++;
	}

	kmer_database.Close();

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	cerr << "Master kmer database processed: " << duration.count() << " seconds" << endl;

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

	// get kmer databases from file
	vector<CharString> kmer_dbs = createFileList(options.kmerDatabasesFilenames);

	// create output array from master kmerdb
	vector<kmer> pa_matrix;
	uint32 kmer_length;
	createOutputArray(options.masterKmerDatabase, length(kmer_dbs), pa_matrix, kmer_length);

	cerr << "There are " << length(kmer_dbs) << " databases to process." << endl;

	int counter = 0;
	for(auto i : kmer_dbs)
	{
		cout << i << "\t";
		cerr << "Processing " << i << " " << (counter+1) << "/"<< length(kmer_dbs) << endl;
		auto start = high_resolution_clock::now();
		checkForPA(options.masterKmerDatabase, i, pa_matrix, counter);
		auto stop = high_resolution_clock::now();
	        auto duration = duration_cast<seconds>(stop - start);
        	cerr << "Completed in : " << duration.count() << " seconds" << endl;		
		counter++;
	}
	cout << endl;

	for(auto i : pa_matrix)
	{
		//string k;
		//decode(i.k, k, kmer_length);
		//cout << k << "\t";
		cout << i.k << "\t" << i.bits << endl;;
	}

	return 0;

}
