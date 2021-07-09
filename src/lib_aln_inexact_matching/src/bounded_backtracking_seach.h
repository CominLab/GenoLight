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

#ifndef _BACKTRACKING_SEARCH_H
#define _BACKTRACKING_SEARCH_H

#include "occ.h"
#include "sa.h"
#include "pileup.h"

#ifndef BWAIDX_T
#define BWAIDX_T

typedef struct
{
	bwt_t *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;

#endif

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

typedef struct
{
	ubyte_t *seq;
	bool* fix_pos;
	uint32_t len;
	uint8_t max_offset;
	def_param* dp;
	uint8_t offset;
} input_query;

typedef struct
{
	uint64_t w;
	int bid;
} bwt_width_t;

typedef struct
{
	uint32_t info; //i
	uint8_t n_mm :8;
	int last_diff_pos;
	uint64_t k, l; // (k,l) is the SA region of [i,n-1]
} entry_t;

typedef struct
{
	uint32_t n_entries_substack, m_entries;
	entry_t *start_substack;
} substack_t;

typedef struct
{
	int n_sub_stacks, best, n_entries;
	substack_t *stacks;
} stack_t;

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

#ifndef AMB_KMER
#define AMB_KMER

typedef struct
{
		uint32_t num_occur;
		uint32_t* positions_to_ref;
		uint8_t* posAlternativeSnp;
} snp_kmer;

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

#ifdef __cplusplus
extern "C"
{
#endif

	input_query* init_bwt_seq(uint8_t* pattern_to_search, bool* fix_pos, const size_t pattern_len, const uint8_t max_offset, uint8_t offset, const def_param* dp);

	uint8_t* create_pattern_to_search(const char* pattern_input, const size_t pattern_len);

	void get_approximate_match(input_query* info_seq, output* out, uint16_t* m_hit_found_for, uint16_t* m_hit_found_rev);

	hit_reassembly* get_exact_match_to_compleate_dictionary(bwaidx_t *idx, ubyte_t* seq, int len);

#ifdef __cplusplus
}
#endif

#endif
