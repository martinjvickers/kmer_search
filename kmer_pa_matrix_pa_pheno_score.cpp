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
#include <set>
#include <seqan/seq_io.h>

using namespace seqan2;
using namespace std;
using namespace std::chrono;

struct ModifyStringOptions 
{
	CharString kmerDatabasesFilenames;
	CharString phenotypeFilename;
	int kmer_size;
	CharString outputFilename;
	bool createPAfile;
	CharString geneFilename;
};

struct kmer
{
	uint64 k;
	vector<uint64> bits;
};

/*
 * Inputs:
 * 	list of accessions 	= unified_list.txt
 * 	kmer size		= 31
 * 	phenotype list		= watkins_wheatblast_phenotypes.txt
 *	list of PA matrix files = matrix_list.txt
 *	population scores	= 
 *
 * */
seqan2::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                                    int argc, char const ** argv) {
  ArgumentParser parser("kmer_pa_matrix");
  addOption(parser, ArgParseOption("k", "kmer-databases",
			           "A file with a list of kmer databases",
                                    ArgParseArgument::STRING, "TEXT", true));
  setRequired(parser, "kmer-databases");

  addOption(parser, ArgParseOption("p", "phenotype-file",
			           "Your phenotype file for each accession. Accession names must match the accession files exactly.",
                                   ArgParseArgument::STRING, "TEXT", true));

  addOption(parser, ArgParseOption("s", "kmer-size",
	                           "At the moment we specify the kmer size, this should match the kmer size of the database",
	                            ArgParseArgument::INTEGER, "INT"));
  setRequired(parser, "kmer-size");

  addOption(parser, ArgParseOption("i", "input-matrix-file",
			           "This is the PA matrix file",
				   ArgParseArgument::STRING, "TEXT", true));
  setRequired(parser, "input-matrix-file");



  addOption(parser, ArgParseOption("g", "gene-of-interest",
                                   "This is a fasta file of your gene.",
                                   ArgParseArgument::STRING, "TEXT", true));

  setShortDescription(parser, "kmer_pa_matrix");
  setVersion(parser, "0.0.1");
  setDate(parser, "April 2024");
  addUsageLine(parser, "-k kmerdb -i pa_matrix.bin \
                        [\\fIOPTIONS\\fP] ");
  addDescription(parser, "");
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  if (res != ArgumentParser::PARSE_OK){
    return res;
  }

  getOptionValue(options.kmerDatabasesFilenames, parser, "kmer-databases");
  getOptionValue(options.phenotypeFilename, parser, "phenotype-file");
  getOptionValue(options.outputFilename, parser, "input-matrix-file");
  getOptionValue(options.geneFilename, parser, "gene-of-interest");
  getOptionValue(options.kmer_size, parser, "kmer-size");

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

int encode(const Dna5String kmer, unsigned long long int& x_ull)
{
	x_ull = 0;
	for (size_t i = 0; i < length(kmer); ++i)
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
	uint64 count = 0; // the row we are editing

	CharString lockExtension = ".lock.";
	CharString lockFilename = inputFilename;
	CharString byte_suffix = to_string(byte_to_edit);
	append(lockFilename, lockExtension);
	append(lockFilename, byte_suffix);
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

		position = position + ((numUint64+1)*8);
		count++;
	}

	infile.close();

	deleteLockFile(lockFilename);

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

int readPhenotypeList(CharString phenotypesFilename,
		      map<CharString, vector<double>> &v,
		      vector<CharString> &phenotypeNames)
{
	string line;
	ifstream file(toCString(phenotypesFilename));
	if (file.is_open())
	{
		int count = 0;
		while (getline(file, line))
		{
			//cout << line << endl;
			std::vector<std::string> tokens;
			std::istringstream iss(line);
			std::string token;
			while(std::getline(iss, token, '\t'))
			{
				tokens.push_back(token);
			}

			// list of phenotype names
			if(count == 0)
			{
				for(int i = 1; i < tokens.size(); i++)
				{
					phenotypeNames.push_back(tokens[i]);
				}
			}
			else // the phenotype scores for each accession
			{	
				vector<double> scores;
				for(int i = 1; i < tokens.size(); i++)
				{
					scores.push_back(stod(tokens[i]));
				}
				v.insert({tokens[0], scores});
			}
			count++;
		}
		file.close();

	}
	return 0;
}

void createKmerFileList(CharString kmers_to_search_for, std::set<uint64> &v)
{
	string line;
	ifstream file(toCString(kmers_to_search_for));
	if (file.is_open())
	{
		while (getline(file, line))
		{
			unsigned long long int encoded;
	                encode(line, encoded);
			v.insert(encoded);
		}
		file.close();
	}
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

uint64_t reverseBits(uint64_t n)
{
	bitset<64> bits(n);
	string s = bits.to_string();
	std::reverse(s.begin(), s.end());
	return bitset<64>(s).to_ulong();
}

int readInPA(vector<kmer> &pa_matrix, CharString inputFilename, vector<CharString> kmer_dbs, std::set<uint64> &to_find, int kmer_size)
{
        std::ifstream infile(toCString(inputFilename), std::ios::binary);
	if (!infile.is_open())
	{
		std::cerr << "Failed to open file: " << inputFilename << std::endl;
		return 1;
	}

	kmer matrix;
	size_t numUint64 = (length(kmer_dbs) + 63) / 64;
	int bits_at_the_end = 64 - ((numUint64 * 64) - length(kmer_dbs));

	// read a kmer
	while (infile.read(reinterpret_cast<char*>(&matrix.k), sizeof(matrix.k)))
	{
		string kmer;
		decode(matrix.k, kmer, kmer_size); 

		if(to_find.find(matrix.k) != to_find.end())
		{
			cout << kmer << "\t" << matrix.k << "\t";
		}

		for (size_t i = 0; i < numUint64; ++i)
		{
			uint64 byte;
			infile.read(reinterpret_cast<char*>(&byte), sizeof(uint64));
			uint64 rev = reverseBits(byte);
			std::bitset<64> x(rev);
			if(to_find.find(matrix.k) != to_find.end())
	                {
				if(i < numUint64-1)
				{
					cout << x;
				}
				else if(i == numUint64-1)
				{
					for(int j = 0; j < bits_at_the_end; j++)
					{
						std::bitset<64> t(byte);
						cout << t[j];
					}
				}
			}
		}

		if(to_find.find(matrix.k) != to_find.end())
		{
			cout << endl;
		}

	}
	infile.close();
	return 0;
}

void kmers_from_gene(CharString filename, std::set<uint64> &kmers_to_search_for, int kmer_size)
{
	CharString id;
	Dna5String seq;
	SeqFileIn seqFileIn(toCString(filename));
	readRecord(id, seq, seqFileIn);

	for(int i = 0; i < length(seq)-kmer_size-1; i++)
	{
		Dna5String inf = infix(seq, i, i+kmer_size);
		unsigned long long int encoded;
                encode(inf, encoded);
		kmers_to_search_for.insert(encoded);
		reverseComplement(inf);
		encode(inf, encoded);
		kmers_to_search_for.insert(encoded);
	}
	close(seqFileIn);
/*
	cout << "Kmers to search for" << endl;
	for(auto i : kmers_to_search_for)
	{
		string kmer;
		decode(i, kmer, kmer_size);
		cout << kmer << "\t" << i << endl;
	}
*/
}

int openAllFiles(vector<ifstream> &fileStreams, int numFiles, vector<CharString> &filenames)
{
	// Open all files
	for (int i = 0; i < numFiles; i++)
	{
		fileStreams[i].open(toCString(filenames[i]), std::ios::binary);
		if (!fileStreams[i])
		{
			std::cerr << "Error opening file: " << filenames[i] << std::endl;
			return 1; // Exit if any file fails to open
		}
	}
	return 0;
	
}

int closeAllFiles(vector<ifstream> &fileStreams, int numFiles, vector<CharString> &filenames)
{
	// Open all file
	for (int i = 0; i < numFiles; i++)
	{
		fileStreams[i].close();
		cout << "Closed\t" << filenames[i] << endl;
	}
	return 0;
}	

int readLineFromAllFilesBuf(vector<ifstream> &fileStreams, int numFiles, vector<kmer> &matrix_buf_vec, int mb)
{
	// 1MB buffer size
	const size_t buffer_size = 1 * 1024 * 1024;
	//const size_t buffer_size = 64;
	
        for(int i = 0; i < numFiles; i++)
	{
		vector<uint64_t> buffer(buffer_size);
		
		size_t read_count = fileStreams[i].read(reinterpret_cast<char*>(buffer.data()), buffer_size * sizeof(uint64_t)).gcount() / sizeof(uint64_t);
		if (read_count == 0) break;  // End of file

		kmer matrix;
		int kmer_counter = 0;

		for(int j = 0; j < read_count; j+=2)
		{
			uint64_t k = buffer[j];
			uint64_t b = buffer[j+1];

			if(i == 0) // populate the matrix_buf_vec if not done before
			{
				kmer result;
				result.k = k;
				vector<uint64_t> vec(numFiles, 0);
				vec[i] = b;
				result.bits = vec;
				
				matrix_buf_vec.push_back(result);
			}
			else
			{
				matrix_buf_vec[kmer_counter].bits[i] = b;
			}

			kmer_counter++;
		}
	}

        return 0;
}

/* process phenotype scores
 *    take the matrix_buf and the names of each accession
 *    I have a kmer vector AGT, 000011010
 *    
 */
int process_phenotype_scores(vector<kmer> &matrix_buf, vector<CharString> &accessionNamesInPA,
			     map<CharString, vector<double>> &phenotypes, int num_phenoms,
			     vector<int> &pheno_to_accession_map)
{
	size_t numUint64 = (accessionNamesInPA.size() + 63) / 64;

	vector<double> vec(num_phenoms, 0.0);

	//loop through each kmer in the buffer
	for(auto k : matrix_buf)
	{
		int num_bits_set = 0;
		fill(vec.begin(), vec.end(), 0.0); // reset values back to zero

		// loop through each accession with a phenotype
		for(auto i : pheno_to_accession_map)
		{
			//cout << accessionNamesInPA[i] << " exists" << endl;
			auto &phenotype_vec = phenotypes[accessionNamesInPA[i]];

			// if so, loop through the phenotyp to sum the score
			//int byte_to_edit = (((i+1) + 63) / 64)-1;
			//int bit_to_edit = i - (byte_to_edit*64);

			int byte_to_edit = i / 64;
			int bit_to_edit = i % 64;

			//cout << accessionNamesInPA[i] << " exists at position " << i << " byte " << byte_to_edit << " " << bit_to_edit;
			// is this bit set?
			if (k.bits[byte_to_edit] & (1ULL << bit_to_edit))
			{
				num_bits_set++;
				//cout << "bit is set ";
				for(int p = 0; p < num_phenoms; p++)
				{
					//vec[p] = vec[p] + (1 * phenotypes[accessionNamesInPA[i]][p]);
					vec[p] = vec[p] + (1 * phenotype_vec[p]);
					//	cout << vec[p] << " ";
				}
			}
				//cout << endl;
		}
		
/*
		string kmer;
		decode(k.k, kmer, 31);
		cout << kmer << "\t";
		for(auto s : vec)
		{
			cout << s << "\t";
			cout << s/num_bits_set << "\t";
		}
		cout << num_bits_set << "\t";
		cout << endl;
*/
	}

	return 0;
}

int work(vector<ifstream> &fileStreams, vector<CharString> &matrixFilenames, map<CharString, vector<double>> &phenotypes, 
	 vector<CharString> &kmer_dbs, vector<CharString> &phenotypeNames, vector<int> &pheno_to_accession_map)
{
	int counter = 0;
	int num_mb = 1000; //AKA 1GB
	vector<kmer> matrix_buf;
	do
	{
		auto start = high_resolution_clock::now();
		matrix_buf.clear();
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);
		cout << "Clearing previous buffer " << duration.count() << endl;
		// test break while working
		//if(counter >= 4)
		//	break;

		start = high_resolution_clock::now();
		readLineFromAllFilesBuf(fileStreams, matrixFilenames.size(), matrix_buf, num_mb);
		stop = high_resolution_clock::now();
		duration = duration_cast<seconds>(stop - start);
		cout << "Read 1 block of 1MB in " << duration.count() << endl;

		// print the first one
		/*
		cout << matrix_buf[0].k << "\t";
		for(auto i : matrix_buf[0].bits)
		{
			cout << i << "\t";
		}
		cout << endl;
		*/

		// now we actually process something
		start = high_resolution_clock::now();
		process_phenotype_scores(matrix_buf, kmer_dbs, phenotypes, phenotypeNames.size(), pheno_to_accession_map);
		stop = high_resolution_clock::now();
		duration = duration_cast<seconds>(stop - start);
		cout << "Time to process that block " << duration.count() << endl;

		counter++;

	} while(matrix_buf.size() > 0);

	return 0;
}

// this function creates and opens the file streams, initiates the work and then closes them.
int process(vector<CharString> &matrixFilenames, vector<CharString> &kmer_dbs,
	    map<CharString, vector<double>> &phenotypes, vector<CharString> &phenotypeNames,
	    vector<int> &pheno_to_accession_map)
{
	// open files
	vector<ifstream> fileStreams;
	fileStreams.resize(matrixFilenames.size());
	openAllFiles(fileStreams, matrixFilenames.size(), matrixFilenames);

	work(fileStreams, matrixFilenames, phenotypes, kmer_dbs, phenotypeNames, pheno_to_accession_map);

	cout << "finished" << endl;

	// close files
	closeAllFiles(fileStreams, matrixFilenames.size(), matrixFilenames);

	return 0;
}

vector<int> createPheno2AM(vector<CharString> accession_names, map<CharString, vector<double>> phenotypes)
{
	vector<int> pheno_to_accession_map;

	for(int a = 0; a < accession_names.size(); a++)
	{	
		if (phenotypes.find(accession_names[a]) != phenotypes.end())
		{
			pheno_to_accession_map.push_back(a);
		}
	}

	return pheno_to_accession_map;
}

// A basic template to get up and running quickly
// ./kmer_pa_matrix_search -k examples/massive_list.txt -i meh.out -l kmer.txt -s 3
//
//
// ./kmer_pa_matrix_pa_pheno_score -k /mnt/platforms/informatics/mvickers/watkins_PA_matrix/unified_list.txt -s 31 -p /mnt/platforms/informatics/mvickers/watkins_PA_matrix/watkins_wheatblast_phenotypes.txt -i /mnt/platforms/informatics/mvickers/watkins_PA_matrix/matrix_list.txt
int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
	if(res != ArgumentParser::PARSE_OK){
		return res == ArgumentParser::PARSE_ERROR;
	}

	// get list of accessions
	vector<CharString> accession_names = createFileList(options.kmerDatabasesFilenames);

	// get phenotypes
	
	map<CharString, vector<double>> phenotypes;
	vector<CharString> phenotypeNames;
	readPhenotypeList(options.phenotypeFilename, phenotypes, ref(phenotypeNames));


	// what we need is a vector of mappings from the phenotype name to the location in kmer_dbs
	vector<int> pheno_to_accession_map = createPheno2AM(accession_names, phenotypes);

	/*
	for(auto i : phenotypeNames)
		cout << i << endl;
	for(auto i : phenotypes)
		cout << i.first << "\t" << i.second[0] << endl;
	*/

	// get list of matrix files
	vector<CharString> matrixFilenames = createFileList(options.outputFilename);

	process(matrixFilenames, accession_names, phenotypes, phenotypeNames, pheno_to_accession_map);

	return 0;

}
