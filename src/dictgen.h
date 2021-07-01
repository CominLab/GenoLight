#ifndef DICTGEN_H
#define DICTGEN_H

#include <fstream>
#include <iostream>
#include "ParseFasta.hpp"
#include "lib_aln_inexact_matching.h"

#define SNP_INFO_MAKE(pos, ref) ((((pos) & 0x1F) << 3) | ((ref) & 0x07))
#define SNP_INFO_POS(snp_info)  (((snp_info) & 0xF8) >> 3)
#define SNP_INFO_REF(snp_info)  ((snp_info) & 0x07)

/* table for storing multiple positions in case of ambiguous k-mers */
#define AUX_TABLE_COLS 10

#define AUX_TABLE_INIT_SIZE 75000000

#define POS_AMBIGUOUS ((uint32_t)(-1))
#define FLAG_UNAMBIGUOUS 0x00
#define FLAG_AMBIGUOUS   0x01
#define SNP_AUX_TABLE_INIT_SIZE 10000000


typedef uint64_t kmer_t; /* k = 32 */
typedef uint8_t snp_info; /* 5 bits: pos in kmer, 3 bits: ref base */

struct aux_table
{
		uint32_t pos_list[AUX_TABLE_COLS];
}__attribute__((packed));


struct kmer_info
{
		kmer_t kmer;
		uint32_t pos;
}__attribute__((packed));

struct snp_kmer_info
{
	kmer_t kmer;
	uint32_t pos;
	snp_info snp;
	uint8_t alt_base;
	uint8_t ref_freq;
	uint8_t alt_freq;
}__attribute__((packed));

struct snp_aux_table_dictgen
{
	kmer_t kmer;
	uint32_t pos_list[AUX_TABLE_COLS];
	snp_info snp_list[AUX_TABLE_COLS];
	uint8_t alt_base[AUX_TABLE_COLS];
	uint8_t ref_freqs[AUX_TABLE_COLS];
	uint8_t alt_freqs[AUX_TABLE_COLS];
}__attribute__((packed));

struct param_creation_dict
{
		uint8_t k;
		std::string ref_filename;
		std::string snp_filename;
		bwaidx_t* index_reassembly;
		std::string path_info_file;
		std::string chrlens_filename;
		std::string path_ref_genome_index;
		std::string path_reassembly_snp;
		std::string snp_dict_filename;
		std::string snpTemp_filename;
};

void make_snp_dict_from_vcf(ParseFasta genome, std::string snp_dict_filename, std::string snp_filename, std::string snpTemp_filename, bool** snp_locations, size_t size);

param_creation_dict vargeno_init(int argc, char* argv[]);

#endif /* DICTGEN_H */

