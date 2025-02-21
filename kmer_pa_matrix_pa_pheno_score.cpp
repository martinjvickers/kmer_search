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
	int threads = 1;
	int blocks = 1;
};

ModifyStringOptions options; // global

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

  addOption(parser, ArgParseOption("t", "threads",
			           "Number of threads to use. Ideally the same as the number of cpu-cores you have available to you.",
				   ArgParseArgument::INTEGER, "INT"));

  addOption(parser, ArgParseOption("b", "blocks",
			           "Number of blocks of 1MB to read into memory for each PA matrix you have",
				   ArgParseArgument::INTEGER, "INT"));

  setShortDescription(parser, "kmer_pa_matrix");
  setVersion(parser, "0.0.1");
  setDate(parser, "Feb 2025");
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
  getOptionValue(options.kmer_size, parser, "kmer-size");
  getOptionValue(options.threads, parser, "threads");
  getOptionValue(options.blocks, parser, "blocks");

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

uint64_t reverseBits(uint64_t n)
{
	bitset<64> bits(n);
	string s = bits.to_string();
	reverse(s.begin(), s.end());
	return bitset<64>(s).to_ulong();
}

int openAllFiles(vector<ifstream> &fileStreams, int numFiles, vector<CharString> &filenames)
{
	// Open all files
	for (int i = 0; i < numFiles; i++)
	{
		fileStreams[i].open(toCString(filenames[i]), std::ios::binary);
		if (!fileStreams[i])
		{
			cerr << "Error opening file: " << filenames[i] << std::endl;
			return 1; 
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
		cerr << "Closed\t" << filenames[i] << endl;
	}
	return 0;
}	

int readLineFromAllFilesBuf(vector<ifstream> &fileStreams, int numFiles, vector<kmer> &matrix_buf_vec, int mb)
{
	// 1MB buffer size
	const size_t buffer_size = 1 * 1024 * 1024;

	// work on one file
        for(int i = 0; i < numFiles; i++)
	{
		int kmer_counter = 0;

		// populate matrix_buf_vec with mb X buffersize of data
		for(int m = 0; m < mb; m++)
		{
			vector<uint64_t> buffer(buffer_size);
		
			size_t read_count = fileStreams[i].read(reinterpret_cast<char*>(buffer.data()), buffer_size * sizeof(uint64_t)).gcount() / sizeof(uint64_t);
			if (read_count == 0) break;  // End of file

			// add everything to memory
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
	}

        return 0;
}

struct ThreadResult
{
	kmer k_value;
	vector<double> vec;
	int num_bits_set;
};

struct ThreadData
{
	vector<kmer> *matrix_buf;
	vector<CharString> *accessionNamesInPA;
	map<CharString, vector<double>> *phenotypes;
	vector<int> *pheno_to_accession_map;
	int num_phenoms;
	size_t start_idx, end_idx;
	vector<ThreadResult> *results; // Store results
	pthread_mutex_t *mutex; // Synchronization
};

void *process_chunk(void *arg)
{
	ThreadData *data = static_cast<ThreadData *>(arg);
	vector<kmer> &matrix_buf = *data->matrix_buf;
	vector<CharString> &accessionNamesInPA = *data->accessionNamesInPA;
	map<CharString, vector<double>> &phenotypes = *data->phenotypes;
	vector<int> &pheno_to_accession_map = *data->pheno_to_accession_map;
	int num_phenoms = data->num_phenoms;

	vector<double> vec(num_phenoms, 0.0);
	vector<ThreadResult> local_results; // Store local results

    	for(size_t idx = data->start_idx; idx < data->end_idx; ++idx) 
	{
        	auto &k = matrix_buf[idx];
        	fill(vec.begin(), vec.end(), 0.0); // Reset vec per k-mer
        	int num_bits_set = 0;

        	for(auto i : pheno_to_accession_map)
		{
			CharString &accession = accessionNamesInPA[i];
			auto &phenotype_vec = phenotypes[accession];

			int byte_to_edit = i / 64;
			int bit_to_edit = i % 64;

			if(k.bits[byte_to_edit] & (1ULL << bit_to_edit))
			{
                		num_bits_set++;
                		for(int p = 0; p < num_phenoms; p++)
				{
                    			vec[p] += phenotype_vec[p];
                		}
            		}
        	}
		//Check if vec contains values > 0
		if (any_of(vec.begin(), vec.end(), [](double v) { return v > 0.01; })) 
		{
			local_results.push_back({k, vec, num_bits_set});
		}
	}

	// Lock and merge local results into the shared results
	pthread_mutex_lock(data->mutex);
	data->results->insert(data->results->end(), local_results.begin(), local_results.end());
	pthread_mutex_unlock(data->mutex);

	pthread_exit(nullptr);
}

int process_phenotype_scores(vector<kmer> &matrix_buf, vector<CharString> &accessionNamesInPA,
		             map<CharString, vector<double>> &phenotypes, int num_phenoms,
			     vector<int> &pheno_to_accession_map, int threads, vector<ThreadResult> &results)
{
	const int NUM_THREADS = options.threads;
	pthread_t mythreads[NUM_THREADS];
	ThreadData thread_data[NUM_THREADS];
	pthread_mutex_t mutex;
	pthread_mutex_init(&mutex, nullptr);

	size_t chunk_size = matrix_buf.size() / NUM_THREADS;

	for (int t = 0; t < NUM_THREADS; ++t)
	{
        	thread_data[t] = {
			&matrix_buf, &accessionNamesInPA, &phenotypes, &pheno_to_accession_map,
            		num_phenoms, t * chunk_size, (t == NUM_THREADS - 1) ? matrix_buf.size() : (t + 1) * chunk_size,
	    		&results, &mutex
        	};
        	pthread_create(&mythreads[t], nullptr, process_chunk, &thread_data[t]);
    	}

    	for (int t = 0; t < NUM_THREADS; ++t)
	{
        	pthread_join(mythreads[t], nullptr);
    	}

    	pthread_mutex_destroy(&mutex);

    	return 0;
}

int printResult(vector<ThreadResult> &results, vector<int> &pheno_to_accession_map)
{
	for(auto r : results)
	{
		// decode and print kmer
		string kmer;
		decode(r.k_value.k, kmer, 31);

		cout << kmer << "\t";

		// print the phenotype scores
		for(auto v : r.vec)
		{
			cout << v << "\t";
		}

		// print the number of bits that were set
		cout << r.num_bits_set<< "\t";

		// experimental:
		// print out the PA bits of ONLY the accession in the
		// phenotype file
                for(auto i : pheno_to_accession_map)
		{
			int byte_to_edit = i / 64;
			int bit_to_edit = i % 64;
		}

		// print the whole PA bits
		for(auto b : r.k_value.bits)
		{
			//get bitset
			uint64 rev = reverseBits(b);
			//reverse it to little endian for readability
			std::bitset<64> x(rev);
			cout << x;
		}

		cout << endl;
	}
	return 0;
}

int work(vector<ifstream> &fileStreams, vector<CharString> &matrixFilenames, map<CharString, vector<double>> &phenotypes, 
	 vector<CharString> &kmer_dbs, vector<CharString> &phenotypeNames, vector<int> &pheno_to_accession_map)
{
	int counter = 0;
	int num_mb = options.blocks; 
	vector<kmer> matrix_buf;
	vector<ThreadResult> results;
	do
	{
		auto start = high_resolution_clock::now();
		matrix_buf.clear();
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);
		cerr << "Clearing previous buffer " << duration.count() << "s" << endl;
		// test break while working
		//if(counter >= 4)
		//	break;

		// read from your files in num_mb blocks of 1MB
		start = high_resolution_clock::now();
		readLineFromAllFilesBuf(fileStreams, matrixFilenames.size(), matrix_buf, num_mb);
		stop = high_resolution_clock::now();
		duration = duration_cast<seconds>(stop - start);
		cerr << "Read " << num_mb << " block(s) of 1MB in " << duration.count() << "s" << endl;


		// now we actually process something
		start = high_resolution_clock::now();
		process_phenotype_scores(matrix_buf, kmer_dbs, phenotypes, phenotypeNames.size(), pheno_to_accession_map, 4, results);
		stop = high_resolution_clock::now();
		duration = duration_cast<seconds>(stop - start);
		cerr << "Time to process that block " << duration.count() << "s " << results.size() << endl;

		printResult(results, pheno_to_accession_map);

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

	cerr << "Finished" << endl;

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
	//ModifyStringOptions options; //I've made options global
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

	// get list of matrix files
	vector<CharString> matrixFilenames = createFileList(options.outputFilename);

	process(matrixFilenames, accession_names, phenotypes, phenotypeNames, pheno_to_accession_map);

	return 0;
}
