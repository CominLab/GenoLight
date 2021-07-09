#include <cstdint>
#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <bitset>
#include <algorithm>
#include <cstdlib>
#include <experimental/filesystem>
#include <sys/stat.h>
#include "utility.h"
#include "ParseFasta.hpp"
#include <assert.h>

/**
 *
 */
void arg_check(int argc, int expected, std::string command)
{
	if (argc - 2 != expected) /* -2 because option params start at argv[2] */
	{
		if (command == "index")
			print_help_building_dictionaries();
		else
			print_help_geno();

		exit (EXIT_FAILURE);
	}
}

/**
 *@path_info_file
 *@len_reassembly_snp
 */
void CreateInfoFile(std::string path_info_file, uint64_t len_reassembly_snp)
{
	std::ofstream info_file_write(path_info_file, std::ios::out);
	info_file_write << "The dictionaries were built with the following parameters:" << std::endl;
	info_file_write << "k : " << 32 << std::endl;
	info_file_write << ":" << std::endl;
	info_file_write << "Len reassembly reference genome modified with the SNPs: " << unsigned(len_reassembly_snp) << std::endl;
	info_file_write.close();
}

/**
 * @path
 */
void create_default_out_folder(std::string path)
{
	//std::string dir = "";
	//std::string delimiter = "/";
	//size_t pos = path.find(delimiter);
	//while (pos != path.npos)
	//{
	//	dir += path.substr(0, pos + delimiter.length());
	//	path.erase(0, pos + delimiter.length());
	//	mkdir(dir.c_str(), S_IRWXU);
	//	pos = path.find(delimiter);
	//}
	mkdir(path.c_str(), S_IRWXU);
}

/**
 * @path
 */
std::string get_name_from_genome_path(std::string path)
{
	std::vector < std::string > s = split_vargeno(path, '/');
	size_t number_element = s.size();
	std::string name_genome = s[number_element - 1];
	return name_genome.substr(0, name_genome.find("."));
}

std::string get_file_from_path(std::string path)
{
	std::vector < std::string > s = split_vargeno(path, '/');
	size_t number_element = s.size();
	return s[number_element - 1];
}

/**
 *
 */
void print_help_valid_comand(void)
{
	std::cerr << "\nPLEASE CHECK the README file to understand how to use GenoLight!\n\n";
	std::cerr << "\n-Valid commands:\n\n";
	std::cerr << "\n(should be used in this order)\n\n";
	std::cerr << "-->create_incomplete_smartSnpDictionary -r [reference genome] -s [snp dictionary] -p [output files prefix]\n" << std::endl;
	std::cerr << "-->reassembly -n [snp dictionary] -p [output files prefix]\n" << std::endl;
	std::cerr << "-->createFMDIndex\n" << std::endl;
	std::cerr << "-->complete_smartSnpDictionary -n [snp dictionary] -p [output files prefix]\n" << std::endl;
	std::cerr << "-->geno -t [number of threads] -r [single end reads] -s [snp dictionary] -p [output files prefix] -o [genotyping output file]\n" << std::endl;
	exit (EXIT_FAILURE);
}

void print_help_incDic(void)
{
	std::cerr << "\nPLEASE CHECK the README file to understand how to use GenoLight!\n\n";
	std::cerr << "\n'create_incomplete_smartSnpDictionary' command needs\n\n";
	std::cerr << "-r [reference genome] -s [snp dictionary] -p [output files prefix]\n" << std::endl;
	exit (EXIT_FAILURE);
}

void print_help_reassembly(void)
{
	std::cerr << "\nPLEASE CHECK the README file to understand how to use GenoLight!\n\n";
	std::cerr << "\n'reassembly' command needs\n\n";
	std::cerr << "-n [snp dictionary] -p [output files prefix]\n" << std::endl;
	exit (EXIT_FAILURE);
}

void print_help_compDic(void)
{
	std::cerr << "\nPLEASE CHECK the README file to understand how to use GenoLight!\n\n";
	std::cerr << "\n'complete_smartSnpDictionary' command needs\n\n";
	std::cerr << "-n [snp dictionary] -p [output files prefix]\n" << std::endl;
	exit (EXIT_FAILURE);
}

void print_help_geno(void)
{
	std::cerr << "\nPLEASE CHECK the README file to understand how to use GenoLight!\n\n";
	std::cerr << "\n'geno' command needs\n\n";
	std::cerr << "-t [number of threads] -r [single end reads] -s [snp dictionary] -p [output files prefix] -o [genotyping output file]\n" << std::endl;
	std::cerr << "\nor\n";
	std::cerr << "-t [number of threads] -rl [Left/.1 paired end reads] -rr [Right/.2 paired end reads] -s [snp dictionary] -p [output files prefix] -o [genotyping output file]\n" << std::endl;
	exit (EXIT_FAILURE);
}

void print_help_completeDictionaty(void)
{
	std::cerr << "MANDATORY:\n\n";
	std::cerr << "....\n\n";

	std::cerr << "\nOPTIONAL:\n";

	std::cerr << "In this version there are no optional parameters for this command.\n\n";
	exit (EXIT_FAILURE);
}

/**
 *
 */
void print_help_building_dictionaries(void)
{
	std::cerr << "MANDATORY:\n\n";

	std::cerr << "-r STRING Path to reference genome.\n";
	std::cerr << "-s STRING Path SNPs database.\n";

	std::cerr << "\nOPTIONAL:\n";

	std::cerr << "-p STRING Dictionary prefix.\n\n";
	exit (EXIT_FAILURE);
}

/**
 * @kmer
 * @kmer_had_n
 */
kmer_t encode_kmer(const char *kmer, bool& kmer_had_n)
{
#define KMER_ADD_BASE(x) (encoded_kmer |= (x))

	kmer_t encoded_kmer = 0UL;
	char *base = (char *) &kmer[0];
	for (int i = 0; i < 32; i++)
	{
		encoded_kmer <<= 2;
		switch (*base++)
		{
			case 'A':
			case 'a':
				KMER_ADD_BASE(0UL);
				break;
			case 'C':
			case 'c':
				KMER_ADD_BASE(1UL);
				break;
			case 'G':
			case 'g':
				KMER_ADD_BASE(2UL);
				break;
			case 'T':
			case 't':
				KMER_ADD_BASE(3UL);
				break;
			case 'N':
			case 'n':
				kmer_had_n = true;
				return 0;
				break;
		}
	}

	kmer_had_n = false;
	return encoded_kmer;

#undef KMER_ADD_BASE
}
/*
 kmer_t encode_kmer(const std::string kmer, uint8_t k_lenght, bool& kmer_had_n)
 {
 #define KMER_ADD_BASE(x) (encoded_kmer |= (x))

 kmer_t encoded_kmer = 0UL;

 for (int i = 0; i < 32; i++)
 {
 encoded_kmer <<= 2;
 switch (kmer[i])
 {
 case 'A':
 case 'a':
 KMER_ADD_BASE(0UL);
 break;
 case 'C':
 case 'c':
 KMER_ADD_BASE(1UL);
 break;
 case 'G':
 case 'g':
 KMER_ADD_BASE(2UL);
 break;
 case 'T':
 case 't':
 KMER_ADD_BASE(3UL);
 break;
 case 'N':
 case 'n':
 kmer_had_n = true;
 return 0;
 default:
 std::abort;
 break;
 }
 }

 kmer_had_n = false;
 return encoded_kmer;

 #undef KMER_ADD_BASE
 }
 */

/**
 * @k
 */
uint64_t calculates_mask(uint8_t k)
{
	uint64_t mask = 3UL;
	uint8_t v = 3;

	for (size_t i = 0; i < k - 1; i++)
		mask = ((mask << 2) | v);
	return mask;
}

uint64_t calculates_mask_for_rc(uint64_t ans, uint8_t k)
{
	uint64_t t = calculates_mask(k);
	t = t << (64 - k * 2);
	ans = ans & t;
	ans = ans >> (64 - k * 2);
	return ans;
}

/**
 * @kmer
 * @next_base
 */
kmer_t shift_kmer(const kmer_t kmer, const char next_base)
{
#define KMER_SHIFT(x) ((kmer << 2) | (x))

	switch (next_base)
	{
		case 'A':
		case 'a':
			return KMER_SHIFT(0UL);
		case 'C':
		case 'c':
			return KMER_SHIFT(1UL);
		case 'G':
		case 'g':
			return KMER_SHIFT(2UL);
		case 'T':
		case 't':
			return KMER_SHIFT(3UL);
			break;
	}
	return 0;

#undef KMER_SHIFT
}

/*
 kmer_t shift_kmer(const kmer_t kmer, const char next_base, uint64_t mask)
 {
 #define KMER_SHIFT(x) (((kmer << 2) | x) & mask)

 switch (next_base)
 {
 case 'A':
 case 'a':
 return KMER_SHIFT(0UL);
 case 'C':
 case 'c':
 return KMER_SHIFT(1UL);
 case 'G':
 case 'g':
 return KMER_SHIFT(2UL);
 case 'T':
 case 't':
 return KMER_SHIFT(3UL);
 default:
 std::abort;
 break;
 }
 return 0;

 #undef KMER_SHIFT
 }
 */

/**
 * @chrlens_filename
 * @genome
 */
void CreateChrlensFile(std::string chrlens_filename, ParseFasta& genome)
{
	std::ofstream chrlens(chrlens_filename);
	//asser

	//Memorize the total number of genome bases
	chrlens << genome.GetNumOfGenomeBases() << std::endl;

	//I memorize the number of chromosomes
	chrlens << genome.Size() << std::endl;

	//For each chromosome, I memorize the name and length
	for (size_t i = 0; i < genome.Size(); i++)
	{
		std::string seq = genome.GetSequence(i)->sequence;
		chrlens << genome.GetSequence(i)->name + " " << seq.size() << std::endl;
	}
	chrlens.flush();
	chrlens.close();
}

/**
 * @chrlens_filename
 */
uint32_t GetNumBasesReferenceGenome(std::string chrlens_filename)
{
	std::string line;
	std::ifstream chrlens_file(chrlens_filename);
	if (!chrlens_file.good())
	{
		std::cerr << "Error create: " << chrlens_filename << std::endl;
		exit (EXIT_FAILURE);
	}
	std::getline(chrlens_file, line);
	uint32_t num_bases = std::stol(line);

	chrlens_file.close();

	return num_bases;
}

/**
 * @text
 * @sep
 */
std::vector<std::string> split_vargeno(const std::string &text, char sep)
{
	std::vector < std::string > tokens;
	std::size_t start = 0, end = 0;

	while ((end = text.find(sep, start)) != std::string::npos)
	{
		tokens.push_back(text.substr(start, end - start));
		start = end + 1;
	}
	tokens.push_back(text.substr(start));
	return tokens;
}

void printint(uint64_t x)
{
	std::bitset < 64 > y(x);
	std::cout << y << std::endl;
}
/*
 * https://thispointer.com/how-to-split-a-string-in-c/
 std::string split implementation by using delimiter as a character.
 */
std::vector<std::string> split_line(std::string strToSplit, char c)
{
	std::stringstream ss(strToSplit);
	std::string item;
	std::vector < std::string > splittedStrings;
	while (std::getline(ss, item, c))
	{
		splittedStrings.push_back(item);
	}
	return splittedStrings;
}

/**
 * @e
 * @k_length
 */
std::string decode_kmer_vargeno(uint64_t e, uint8_t k_length)
{
	std::string kmer_decod = "";
	int i;
	for (i = 0; i < k_length; ++i)
	{
		switch (e & 0x0000000000000003)
		{
			case 0:
				kmer_decod = 'A' + kmer_decod;
				break;

			case 1:
				kmer_decod = 'C' + kmer_decod;
				break;

			case 2:
				kmer_decod = 'G' + kmer_decod;
				break;

			case 3:
				kmer_decod = 'T' + kmer_decod;
				break;
		}
		e >>= 2;
	}
	return kmer_decod;
}

/**
 * @kmer
 * @base
 */
unsigned kmer_get_base(const kmer_t kmer, unsigned base)
{
	unsigned s = 2 * base;
	return ((kmer & (0x3UL << s)) >> s);
}

/*

 kmer_t rev_compl(const kmer_t orig, const uint8_t k)
 {
 static uint16_t table[0x10000] = {};
 static bool init = false;
 if (!init)
 {
 uint16_t x;
 uint16_t y;
 for (unsigned i = 0; i < 0x10000; i++)
 {
 x = 0;
 y = i;
 for (unsigned int j = 0; j < 8; j++)
 {
 x <<= 2;
 switch (y & 3)
 {
 case 0:
 x += 3;
 break;
 case 1:
 x += 2;
 break;
 case 2:
 x += 1;
 break;
 case 3:
 break;
 }
 y >>= 2;
 }
 table[i] = x;
 }
 init = true;
 }

 kmer_t ans;
 uint16_t *c_orig = (uint16_t*) (&orig);
 uint16_t *c_ans = (uint16_t*) (&ans);

 for (unsigned i = 0; i < 4; i++)
 {
 c_ans[i] = table[c_orig[3 - i]];
 }

 return calculates_mask_for_rc(ans, k);
 }
 */

/*
 * @orig
 */
kmer_t rev_compl(const kmer_t orig)
{
	static uint16_t table[0x10000] = {};
	static bool init = false;
	if (!init)
	{
		uint16_t x;
		uint16_t y;
		for (unsigned i = 0; i < 0x10000; i++)
		{
			x = 0;
			y = i;
			for (unsigned int j = 0; j < 8; j++)
			{
				x <<= 2;
				switch (y & 3)
				{
					case 0:
						x += 3;
						break;
					case 1:
						x += 2;
						break;
					case 2:
						x += 1;
						break;
					case 3:
						break;
				}
				y >>= 2;
			}
			table[i] = x;
		}
		init = true;
	}

	kmer_t ans;
	uint16_t *c_orig = (uint16_t*) (&orig);
	uint16_t *c_ans = (uint16_t*) (&ans);

	for (unsigned i = 0; i < 4; i++)
	{
		c_ans[i] = table[c_orig[3 - i]];
	}

	return ans;
}

/**
 * @base
 */
uint64_t encode_base(const char base)
{
	switch (base)
	{
		case 'A':
		case 'a':
			return BASE_A;
		case 'C':
		case 'c':
			return BASE_C;
		case 'G':
		case 'g':
			return BASE_G;
		case 'T':
		case 't':
			return BASE_T;
		case 'N':
		case 'n':
			return BASE_N;
		default:
			return BASE_X;
	}
}

/**
 * @line
 */
static int parseLineForMemory(char* line)
{
	int i = strlen(line);
	while (*line < '0' || *line > '9')
		line++;
	line[i - 3] = '\0';
	i = atoi(line);
	return i;
}

/**
 *
 */
int getPeakVirtualMemoryUsed()
{
	//Note: this value is in KB!
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];

	while (fgets(line, 128, file) != NULL)
	{
		if (strncmp(line, "VmPeak:", 7) == 0)
		{
			result = parseLineForMemory(line);
			break;
		}
	}
	fclose(file);
	return result;
}


void serialize_uint64(FILE *out, const uint64_t x)
{
	/* assumes little-endian (or, at least, that the binary
	   dictionary files will be generated by the same machine
	   as the one on which they are used) */
	assert(fwrite(&x, sizeof(uint64_t), 1, out));
}

void serialize_uint32(FILE *out, const uint32_t x)
{
	/* assumes little-endian (or, at least, that the binary
	   dictionary files will be generated by the same machine
	   as the one on which they are used) */
	assert(fwrite(&x, sizeof(uint32_t), 1, out));
}

void serialize_uint8(FILE *out, const uint8_t x)
{
	assert(fwrite(&x, sizeof(uint8_t), 1, out));
}



/**
 * @in
 */
uint64_t read_uint64_t(std::ifstream &in)
{
	uint64_t x;
	in.read(reinterpret_cast<char *>(&x), sizeof(uint64_t));
	return x;
}

/**
 * @in
 */
uint32_t read_uint32_t(std::ifstream &in)
{
	uint32_t x;
	in.read(reinterpret_cast<char *>(&x), sizeof(uint32_t));
	return x;
}

/**
 * @in
 */
uint8_t read_uint8_t(std::ifstream &in)
{
	uint8_t x;
	in.read(reinterpret_cast<char *>(&x), sizeof(uint8_t));
	return x;
}

/**
 * @out
 * @value
 */
void write_uint64_t(std::ostream &out, uint64_t value)
{
	out.write(reinterpret_cast<char *>(&value), sizeof(uint64_t));
}

/**
 * @out
 * @value
 */
void write_uint32_t(std::ostream &out, uint32_t value)
{
	out.write(reinterpret_cast<char *>(&value), sizeof(uint32_t));
}

/**
 * @out
 * @value
 */
void write_uint8_t(std::ostream &out, uint8_t value)
{
	out.write(reinterpret_cast<char *>(&value), sizeof(uint8_t));
}

/**
 *Method that produces a timestamp
 */
std::string debugProgress()
{
	time_t now;
	struct tm nowLocal;
	now = time(NULL);
	nowLocal = *localtime(&now);

	return std::to_string(nowLocal.tm_hour) + ":" + std::to_string(nowLocal.tm_min) + ":" + std::to_string(nowLocal.tm_sec);
}

/**
 * @value
 */
std::string uint64_to_string(uint64_t value)
{
	std::ostringstream os;
	os << value;
	return os.str();
}

std::string returnExtensionFileSnp(std::string path)
{
	std::vector < std::string > columns = split_vargeno(path, '.');
	std::string extension = columns[columns.size() - 1];
	return extension;
}
