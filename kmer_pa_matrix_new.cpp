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
#include <boost/dynamic_bitset.hpp>

using namespace seqan2;
using namespace std;
using namespace std::chrono;

struct ModifyStringOptions 
{
	CharString kmerDatabasesFilenames;
	CharString masterKmerDatabase;
	int specific_kmer;
	CharString outputFilename;
};

struct kmer
{
	uint64 k;
	vector<uint64> bits;
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
  addOption(parser, ArgParseOption("n", "specific-kmer", 
			           "If you wish to just run PA on a specific kmer database, give that number here.",
			           ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("o", "output-file",
			           "This is the PA matrix output file",
				   ArgParseArgument::STRING, "TEXT", true));
  

  setShortDescription(parser, "kmer_pa_matrix");
  setVersion(parser, "0.0.1");
  setDate(parser, "April 2024");
  addUsageLine(parser, "-k kmerdb \
                        [\\fIOPTIONS\\fP] ");
  addDescription(parser, "");
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  if (res != ArgumentParser::PARSE_OK){
    return res;
  }

  getOptionValue(options.kmerDatabasesFilenames, parser, "kmer-databases");
  getOptionValue(options.masterKmerDatabase, parser, "master-database");
  getOptionValue(options.specific_kmer, parser, "specific-kmer");
  getOptionValue(options.outputFilename, parser, "output-file");

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

int checkForPA(CharString currentDBtoCheck, vector<kmer> &v, int n)
{
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
	CKmerAPI kmer_object(_kmer_query_length);

	uint64 integer_to_work_on = ((n+1) + 63) / 64;
	uint64 integer_before = integer_to_work_on-1;
	uint64 num_bits = n - (integer_before*64);

        for(uint64 c = 0; c < length(v); c++)
	{
		string k;
		decode(v[c].k, k, _kmer_query_length);
		kmer_object.from_string(k);
		string meh;
		kmer_object.to_string(meh);
		if(kmer_query_database.IsKmer(kmer_object) == 1)
		{
			v[c].bits[integer_to_work_on-1] ^= (static_cast<uint64_t>(1) << num_bits);
		}
	}

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
	vector<bool> x(n, false);
	size_t bytes_used = x.size() / 8 + (x.size() % 8 != 0);
	cout << "bytes used to represent a kmer: " << bytes_used << endl;
	size_t totalram = ((bytes_used+8)*_total_kmers) / (1024 * 1024 * 1024);
	cout << "estimated RAM " << totalram << " GB" << endl;

	vector<uint64> bits;
	size_t numUint64 = (n + 63) / 64;
	for(int i = 0; i < numUint64; i++)
	{
		uint64 b = 0;
		bits.push_back(b);
	}

	while(kmer_database.ReadNextKmer(kmer_object, counter))
	{
		string str;
		kmer_object.to_string(str);
		unsigned long long int encoded;
		encode(str, encoded);
		v[c].k = encoded;
		v[c].bits = bits;
		c++;
	}

	kmer_database.Close();

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	cerr << "Master kmer database processed: " << duration.count() << " seconds" << endl;

	return 0;
}

int writeOutPA(vector<kmer> &pa_matrix, CharString outputFilename)
{
        // write out
	std::ofstream outfile(toCString(outputFilename), std::ios::binary);
	if (!outfile.is_open())
	{
		std::cerr << "Failed to open file: " << outputFilename << std::endl;
		return 1;
	}

	for(auto i : pa_matrix)
	{
		outfile.write(reinterpret_cast<const char*>(&i.k), sizeof(i.k));
		for (auto bit : i.bits)
		{
			outfile.write(reinterpret_cast<const char*>(&bit), sizeof(bit));
		}
	}
	outfile.close();

	return 0;
}

int readInPA(vector<kmer> &pa_matrix, CharString inputFilename, vector<CharString> kmer_dbs)
{
        std::ifstream infile(toCString(inputFilename), std::ios::binary);
	if (!infile.is_open())
	{
		std::cerr << "Failed to open file: " << inputFilename << std::endl;
		return 1;
	}

	kmer matrix;
	size_t numUint64 = (length(kmer_dbs) + 63) / 64;

	while (infile.read(reinterpret_cast<char*>(&matrix.k), sizeof(matrix.k)))
	{
		cout << matrix.k << "\t";
		for (size_t i = 0; i < numUint64; ++i)
		{
			uint64 bit;
			infile.read(reinterpret_cast<char*>(&bit), sizeof(uint64));
			cout << bit << "\t";
		}
		cout << endl;
	}
	infile.close();
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
		if(counter == options.specific_kmer)
		{
			cout << i << "\t";
			cerr << "Processing " << i << " " << (counter+1) << "/"<< length(kmer_dbs) << endl;
			auto start = high_resolution_clock::now();
			checkForPA(i, pa_matrix, counter);
			auto stop = high_resolution_clock::now();
		        auto duration = duration_cast<seconds>(stop - start);
	        	cerr << "Completed in : " << duration.count() << " seconds" << endl;
		}
		counter++;
	}
	cout << endl;


	writeOutPA(pa_matrix, options.outputFilename);

	//readInPA(pa_matrix, options.outputFilename, kmer_dbs);

	return 0;

}
