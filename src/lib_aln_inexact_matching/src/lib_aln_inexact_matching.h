/* The MIT License

 Copyright (c) 2019 Mattia Marcolin.

 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

/*
 * This header provides public interfaces for using library from codes.
 */

#ifndef LIB_BWA_BACKTRACKING_H_
#define LIB_BWA_BACKTRACKING_H_

//TYPE INDEX
#ifndef TYPE_INDEX_
#define TYPE_INDEX_

#define BWTALGO_AUTO  0
#define BWTALGO_RB2   1
#define BWTALGO_BWTSW 2
#define BWTALGO_IS    3

#endif

#ifndef BNTANN1_T
#define BNTANN1_T

typedef struct
{
		int64_t offset;
		int32_t len;
		int32_t n_ambs;
		uint32_t gi;
		int32_t is_alt;
		char *name, *anno;
} bntann1_t;

#endif

#ifndef BNTAMB1_T
#define BNTAMB1_T

typedef struct
{
		int64_t offset;
		int32_t len;
		char amb;
} bntamb1_t;

#endif

#ifndef BNTSEQ_T
#define BNTSEQ_T

typedef struct
{
		int64_t l_pac;
		int32_t n_seqs;
		bntann1_t *anns; // n_seqs elements
		int32_t n_holes;
		bntamb1_t *ambs; // n_holes elements
		FILE *fp_pac;
} bntseq_t;

#endif

#ifndef BWT_T
#define BWT_T

typedef struct
{
		uint64_t primary; // S^{-1}(0), or the primary index of BWT
		uint64_t L2[5]; // C(), cumulative count
		uint64_t seq_len; // sequence length
		uint64_t bwt_size; // size of bwt, about seq_len/4
		uint32_t *bwt; // BWT
		uint32_t cnt_table[256]; // occurance array, separated to two parts
		// suffix array
		int sa_intv;
		uint64_t n_sa;
		uint64_t *sa;
} bwt_t;

#endif

#ifndef BWAIDX_T
#define BWAIDX_T

typedef struct
{
		bwt_t *bwt; // FM-index
		bntseq_t *bns; // information on the reference sequences
		uint8_t *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;

#endif

#ifndef SEARCH_RESULT
#define SEARCH_RESULT

typedef struct
{
		uint64_t hit;
		uint8_t* altSnpPosInsideKmer;
		uint32_t* positions_to_ref;
		uint32_t num_occur;
		uint32_t* different_positions;
		uint8_t n_mismatches;
		uint8_t offset;
} search_result;

#endif

#ifndef OUTPUT
#define OUTPUT

typedef struct
{
		uint32_t n_for;
		search_result** res_for;

		uint32_t n_rev;
		search_result** res_rev;
} output;

#endif

#ifndef HIT_REASSEMBLY
#define HIT_REASSEMBLY

typedef struct
{
		uint64_t* positions_to_ref;
		bool isRev;

} hit_reassembly;

#endif

#include <stdbool.h>
#include "pileup.h"

#ifndef DEFAULT_PARAM
#define DEFAULT_PARAM

typedef struct
{
		bwaidx_t *idx_reference;
		bwaidx_t *idx_reassembly;
		SNPSmartDictionary* snpSD;
		char quality_threshold;
		uint8_t max_diff;
} def_param;

#endif

#ifdef __cplusplus
extern "C"
{
#endif

	/**
	 *Builds the FMD-index for the reference genome.
	 *
	 *@param path_genome: Path where is locate database sequences in the FASTA format
	 *@param prefix: Prefix of the output database
	 *@param algo_type: Index construction algorithm
	 */
	void lib_aln_index(const char* path_genome, const char* prefix, int algo_type);

	/**
	 *Method that load index in memory.
	 *
	 *@param path_genome:  Path where is locate database sequences in the FASTA format
	 */
	bwaidx_t* lib_aln_idx_load(const char *path_genome);

	/**
	 *Method that free memory by index.
	 *
	 *@param idx: FMD-Index
	 */
	void lib_aln_idx_destroy(bwaidx_t *idx);

	/**
	 *Method that solve approximate pattern matching problem.
	 *
	 * This function dynamically allocates memory. You need to free the memory by
	 * libbwa_bt_sr_destroy().
	 *
	 *@param idx: FMD-Index
	 *@param pattern_input: Input string
	 *@param max_mismatches: Max number of mismatch between hit and reference
	 *@param num_hits: Number of admissible hit found
	 *@param type_search: Defines the type of search. Permissible value are ARBITRARY_HIT and ALL_HITS
	 *@param type_output: Defines the type of output: Permissible value are ALLOW_REV_COMP and NO_REV_COMP
	 */
	output* lib_aln_bound_backtracking(const def_param* dp, const char* pattern_input, const char* fix_pos);


	def_param* lib_aln_set_default_parameters(const bwaidx_t* idx_reference, const bwaidx_t* idx_reassembly, const SNPSmartDictionary* snpSD,
											  const char quality_threshold, const uint8_t max_mismatches);

	hit_reassembly* lib_aln_bound_backtracking_record(const bwaidx_t* idx, const char* kmer);
	/**
	 *Free memory allocate for store the results of the search
	 *
	 *@param result: Output of method lib_aln_bound_backtracking
	 *@param num_hit: Number of hit inside result
	 */
	void lib_aln_sr_destroy(output* result);

	void lib_aln_hit_reassembly_destroy(hit_reassembly* hit_reassembly);

#ifdef __cplusplus
}
#endif

#endif
