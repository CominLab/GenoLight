#ifndef UTILITY_H
#define UTILITY_H

#include <cstdint>
//PER LA GESTIONE DEL PARSER
#include <vector>
#include <limits>
#include <iostream>
#include <fstream>
#include "pileup.h"
#include "lib_aln_inexact_matching.h"
#include "ParseFasta.hpp"
#include "dictgen.h"
#include <iostream>
#include <fstream>

#define POW_2_32 (1UL << 32)
#define UNUSED(x) (void)(x)

//Non ancora usate
#define SET_FLAG_REV(flag) (flag | 0x02)
#define GET_UNICITY(flag) (flag & 0x01)
#define GET_NUM_POS(flag) ((flag & 0x7c) >> 2)

#define BASE_A 0x0
#define BASE_C 0x1
#define BASE_G 0x2
#define BASE_T 0x3
#define BASE_N 0x4
#define BASE_X 0x7
#define CONST_QUALITY_SCORE 33

struct param_creation_geno
{
		std::size_t quality_threshold;
		std::size_t n_threads; //number of threads
		std::size_t chunk_size;
		std::size_t max_offset;
		std::size_t inner_distance;
		std::size_t read_len;
		std::size_t k;

		bwaidx_t* indexReferenceGenome;
		bwaidx_t* indexSNPsReassembly;

		std::string chrlens_filename;
		std::string snp_dict_filename;

		std::string path_snp;
		std::string path_fastq;

		std::string path_fastq_left;
		std::string path_fastq_right;

		std::string path_output;
		std::string path_temp;
		std::string prefix;
		std::string name;
		std::string reference;

		PileupTable ptable;

		//struct packed_pileup_entry *pileup_table;
		uint32_t pileup_size;
};

typedef std::pair<std::size_t, std::size_t> chunk_view;

typedef uint8_t k_t;
const k_t k_max = std::numeric_limits < k_t > ::max();

void arg_check(int argc, int expected, std::string command);

void standard_help(void);

void print_help_valid_comand(void);

void print_help_incDic(void);

void print_help_building_dictionaries(void);

void print_help_compDic(void);

void print_help_completeDictionaty(void);

void print_help_geno(void);

void print_help_reassembly(void);

void printint(uint64_t x);

std::vector<std::string> split_vargeno(const std::string &text, char sep);
std::vector<std::string> split_vargeno(const std::string &text, char sep);
std::vector<std::string> vcf_split_info(const std::string text);
std::string returnExtensionFileSnp(std::string extension);

param_creation_dict* vargeno_init_compleat_dictionary(int argc, char* argv[]);
void CreateChrlensFile(std::string filename, ParseFasta& genome);
uint32_t GetNumBasesReferenceGenome(std::string chrlens_filename);

bool ExtractFrequency(std::string info, float& f1, float& f2);

uint64_t calculates_mask(uint8_t k);
//kmer_t encode_kmer(const char* kmer, bool& kmer_had_n);

//kmer_t encode_kmer(const std::string kmer, uint8_t k_lenght, bool& kmer_had_n);
kmer_t encode_kmer(const char *kmer, bool& kmer_had_n);

kmer_t shift_kmer(const kmer_t kmer, const char next_base);
unsigned kmer_get_base(const kmer_t kmer, unsigned base);
uint64_t encode_base(const char base);
kmer_t rev_compl(const kmer_t orig);

std::string get_name_from_genome_path(std::string path);
std::string get_file_from_path(std::string path);
void create_default_out_folder(std::string path);
void CreateInfoFile(std::string path_info_file, uint64_t len_reassembly_snp);
std::string decode_kmer_vargeno(uint64_t e,uint8_t k_length);

int getPeakVirtualMemoryUsed(void);
std::string uint64_to_string( uint64_t value );

void write_uint64_t(std::ostream &out, uint64_t value);
void write_uint32_t(std::ostream &out, uint32_t value);
void write_uint8_t(std::ostream &out, uint8_t value);

std::string debugProgress();

uint64_t read_uint64_t(std::ifstream &in);
uint32_t read_uint32_t(std::ifstream &in);
uint8_t read_uint8_t(std::ifstream &in);


void serialize_uint64(FILE *out, const uint64_t x);

void serialize_uint32(FILE *out, const uint32_t x);


void serialize_uint8(FILE *out, const uint8_t x);

#endif

