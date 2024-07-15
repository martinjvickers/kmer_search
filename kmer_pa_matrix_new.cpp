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
#include <fstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string>

using namespace seqan2;
using namespace std;
using namespace std::chrono;

struct ModifyStringOptions 
{
	CharString kmerDatabasesFilenames;
	CharString masterKmerDatabase;
	int specific_kmer;
	CharString outputFilename;
	bool createPAfile;
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
  addOption(parser, ArgParseOption("pa", "create-pa",
			           "This creates the PA matrix file, a blank one ready for first use."));

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
  options.createPAfile = isSet(parser, "create-pa");

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

int getListOfKmers(CharString masterDatabase, vector<uint64> &v)
{
	KMC::Runner runner;
	CKMCFile kmer_database;
	if(!kmer_database.OpenForListing(toCString(masterDatabase)))
	{
		cout <<"Error opening kmer database file." << endl;
		return 1;
	}
	// read existing kmer database, which gives us our kmer length
	uint32 _kmer_length, _mode, _counter_size, _lut_prefix_length, _sig_len, _min_count;
	uint64 _max_count;
	uint64 _total_kmers;
	kmer_database.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _sig_len, _min_count, _max_count, _total_kmers);
	CKmerAPI kmer_object(_kmer_length);
	uint32 counter;
	uint64 c;

	v.resize(_total_kmers);

	while(kmer_database.ReadNextKmer(kmer_object, counter))
	{
		unsigned long long int encoded;
		string str;
		kmer_object.to_string(str);
		encode(str, encoded);
		v[c] = encoded;
	}
	
	return 0;
}

int checkForPA(CharString masterKmers, CharString currentKmers, CharString inputFilename, vector<CharString> kmer_dbs, int position_to_edit)
{
	// this vector contains the true/false of each of our master kmers found in the current kmer list
	vector<bool> v;
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
	CKmerAPI kmer_object(_kmer_length);

        // this is the query database for random access searching
	CKMCFile kmer_query_database;
	if(!kmer_query_database.OpenForRA(toCString(currentKmers)))
	{
		cout <<"Error opening kmer database file " << currentKmers << endl;
		return 1;
	}
	uint32 _kmer_query_length, _query_mode, _query_counter_size, _query_lut_prefix_length, _query_sig_len, _query_min_count;
	uint64 _query_max_count, _query_total_kmers;
	kmer_query_database.Info(_kmer_query_length, _query_mode, _query_counter_size, _query_lut_prefix_length, _query_sig_len, _query_min_count, _query_max_count, _query_total_kmers);
	CKmerAPI kmer_query_object(_kmer_query_length);
	uint32 counter;

	// for each master kmer, we append our kmers to the vector
	while(kmer_database.ReadNextKmer(kmer_object, counter))
	{
		unsigned long long int encoded;
		string str;
		kmer_object.to_string(str);
		encode(str, encoded);
		if(kmer_query_database.IsKmer(kmer_object) == 1)
		{
			v.push_back(1);
		}
		else
		{
			v.push_back(0);
		}
	}

	kmer_query_database.Close();
	kmer_database.Close();

	// read in file and edit it
	std::fstream infile(toCString(inputFilename), std::ios::in | std::ios::out | std::ios::binary);
        //if (!infile.is_open())
	if (!infile)
        {
                std::cerr << "Failed to open file: " << inputFilename << std::endl;
                return 1;
        }		

	//this calculates our vector which contains the number of bytes we need to cover our number of accessions
	//this might not be needed. I will give it some thought
	size_t numUint64 = (length(kmer_dbs) + 63) / 64;
	vector<uint64> bits;
	for(int i = 0; i < numUint64; i++)
	{
		uint64 b = 0;
		bits.push_back(b);
	}
	int byte_to_edit = (((position_to_edit+1) + 63) / 64)-1;
	int bit_to_edit = position_to_edit - (byte_to_edit*64);

	// remember that the first byte contains the encoded kmer
	std::streampos position = (0+8) + (byte_to_edit * 8);
	uint64 count = 0;

	while(count < length(v))
	{
		if(v[count] == 1)
		{
			uint64 byte;
			infile.seekg(position);
                	infile.read(reinterpret_cast<char*>(&byte), sizeof(uint64));
		        infile.seekp(position);
			byte = (((uint64_t)1) << bit_to_edit) | byte;
			infile.write(reinterpret_cast<char*>(&byte), sizeof(byte));
		}

		position = (position) + ((numUint64+1)*8);
		count++;
	}

	infile.close();

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

int createOutputArray(CharString masterKmers, int n, vector<kmer> &v, uint32 &kmer_length, CharString outputFilename)
{
        // write out
	std::ofstream outfile(toCString(outputFilename), std::ios::binary);
	if (!outfile.is_open())
	{
		std::cerr << "Failed to open file: " << outputFilename << std::endl;
		return 1;
	}
	
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

	v.resize(_total_kmers);
	uint64 c; // counter

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
		outfile.write(reinterpret_cast<const char*>(&encoded), sizeof(encoded));
		for (auto bit : bits)
		{
			outfile.write(reinterpret_cast<const char*>(&bit), sizeof(bit));
		}
		c++;
	}

	kmer_database.Close();
	outfile.close();

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

// Function to create a lock file with hostname and process ID
bool createLockFile(CharString outputFilename)
{
	int fd = open(toCString(outputFilename), O_CREAT | O_EXCL | O_WRONLY, 0666);
	if (fd == -1)
	{
		return false; // Failed to create lock file, likely because it already exists
	}

	char hostname[256];
        gethostname(hostname, sizeof(hostname));
	pid_t pid = getpid();

	std::string content = "Hostname: " + std::string(hostname) + "\nProcess ID: " + std::to_string(pid) + "\n";
	ssize_t bytesWritten = ::write(fd, content.c_str(), content.size());  // Using ::write to specify from unistd.h

	if(bytesWritten == -1)
	{
		cout << "Error writing to file" << endl;
		close(fd);
		return false;
	}

	if(fsync(fd) == 1)
	{
		cout << "Error syncing to disk" << endl;
		close(fd);
		return false;
	}

	close(fd);
	return true;
}

// Function to delete the lock file
void deleteLockFile(CharString outputFilename)
{
	if (remove(toCString(outputFilename)) != 0) 
	{
		std::cerr << "Error deleting lock file." << std::endl;
	}
	else
	{
		std::cout << "Lock file deleted successfully." << std::endl;
	}
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

	if(options.createPAfile == true)
	{
		createOutputArray(options.masterKmerDatabase, length(kmer_dbs), pa_matrix, kmer_length, options.outputFilename);
	}
	else
	{
		CharString lockExtension = ".lock";
		CharString lockFilename = options.outputFilename;
		append(lockFilename, lockExtension);
		cout << "Lockfile " << lockFilename << endl;
		while(true)
		{
			if(createLockFile(lockFilename))
			{
				break;
			}
			else
			{
				// i guess this could be a random number
				cout << "File locked, waiting 60 seconds" << endl;
				sleep(60);
			}
		}

		cerr << "There are " << length(kmer_dbs) << " databases to process." << endl;
		cout << options.specific_kmer << "\t";
		cerr << "Processing " << options.specific_kmer << " " << (options.specific_kmer+1) << "/"<< length(kmer_dbs) << endl;
		auto start = high_resolution_clock::now();
		checkForPA(options.masterKmerDatabase, kmer_dbs[options.specific_kmer], options.outputFilename, kmer_dbs, options.specific_kmer);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);
		cerr << "Completed in : " << duration.count() << " seconds" << endl;

		// remove lock file
		deleteLockFile(lockFilename);
	}

	//readInPA(pa_matrix, options.outputFilename, kmer_dbs);

	return 0;

}
