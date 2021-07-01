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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "pileup.h"
#include "lib_aln_inexact_matching.h"
#include "bounded_backtracking_seach.h"
#include "fmdindex_load.h"

//Necessary only for GPLv3 version
int fmd_idx_build(const char *fa, const char *prefix, int algo_type, int block_size);

//Derived from bwa_idx_load_from_disk in bwa.c
/**
 *
 */
bwaidx_t* lib_aln_idx_load(const char *path_genome)
{
	bwaidx_t *idx;
	char *prefix;
	prefix = bwa_idx_infer_prefix(path_genome);
	if (prefix == 0)
		return 0;

	idx = calloc(1, sizeof(bwaidx_t));
	idx->bwt = bwa_idx_load_bwt(path_genome);
	if (idx->bwt == 0)
		return 0;
	int i, c;

	idx->bns = bns_restore(prefix);
	if (idx->bns == 0)
		return 0;
	for (i = c = 0; i < idx->bns->n_seqs; ++i)
		if (idx->bns->anns[i].is_alt)
			++c;

	idx->pac = calloc(idx->bns->l_pac / 4 + 1, 1);
	err_fread_noeof(idx->pac, 1, idx->bns->l_pac / 4 + 1, idx->bns->fp_pac); // concatenated 2-bit encoded sequence
	err_fclose(idx->bns->fp_pac);
	idx->bns->fp_pac = 0;

	free(prefix);
	return idx;
}

//Derived from bwa_idx_destroy in bwa.c
/**
 *
 */
void lib_aln_idx_destroy(bwaidx_t *idx)
{
	if (idx == 0)
		return;

	bwt_destroy(idx->bwt);
	//free(idx->bwt);
	free(idx->pac);
	free(idx->bns->anns);
	free(idx->bns->ambs);
	free(idx->bns);
	free(idx);
}

/**
 * This method is used during the phase in which the snp dictionary is completed.
 */
hit_reassembly* lib_aln_bound_backtracking_record(const bwaidx_t *idx_reassembly, const char* kmer)
{
	//Size of pattern input
	size_t pattern_len = strlen(kmer);

	//fprintf(stderr, "\nkmer_to_search: %s", kmer_to_search);

	ubyte_t* pattern_to_search = create_pattern_to_search(kmer, pattern_len);

	//Core method that solves the exact pattern matching problem.
	/*
	 * Since the searched k-mer must be present within the reassembly only once, direct or reverse complement, hit_reassembly is a structure that contains only
	 * the position within the reassembly in which this k-mer was found and if it was found direct or its complementary revere was found
	 */
	hit_reassembly* hit_reassembly = get_exact_match_to_compleate_dictionary(idx_reassembly, pattern_to_search, pattern_len);

	free(pattern_to_search);

	return hit_reassembly;
}

/**
 * This method is invoked for all reads contained in the fastq file.
 */
output* lib_aln_bound_backtracking(const def_param* dp, const char* read_seq, const char* qual)
{
	//Size of pattern input
	size_t pattern_len = strlen(read_seq);

	//From the record I extract the base sequence of the corresponding read and the memorize within this array
	char* kmer_to_search = (char*) malloc(33 * sizeof(char));
	kmer_to_search[32] = '\0';

	bool* fix_position = (bool*) malloc(sizeof(bool) * 32);

	//Calculation of the maximum number of k-mer extractable from the reads
	size_t num_kmer_extractable = pattern_len / 32;
	uint8_t max_offset = pattern_len - (pattern_len % 32) - 32;

	uint16_t m_hit_found_for = 8;
	uint16_t m_hit_found_rev = 8;

	output* out = (output*) malloc(sizeof(output));

	out->n_for = 0;
	out->n_rev = 0;

	out->res_for = (search_result**) malloc(m_hit_found_for * sizeof(search_result*));
	for (int i = 0; i < m_hit_found_for; i++)
		*(out->res_for + i) = (search_result*) malloc(sizeof(search_result));

	out->res_rev = (search_result**) malloc(m_hit_found_rev * sizeof(search_result*));
	for (int i = 0; i < m_hit_found_rev; i++)
		*(out->res_rev + i) = (search_result*) malloc(sizeof(search_result));

	uint8_t* pattern_to_search;
	input_query* seq;

	bool invalidInput = false;

	for (int i = 0; i < num_kmer_extractable; i++)
	{
		uint8_t offset = i * 32;

		for (int m = 0; m < 32; m++)
		{

			if (read_seq[m + offset] != 'N')
			{
				kmer_to_search[m] = read_seq[m + offset];
			}
			else
			{
				invalidInput = true;
				continue;
			}

			kmer_to_search[m] = read_seq[m + offset];

			if (qual[m] > dp->quality_threshold)
				fix_position[32 - m - 1] = true; //TRUE
			else
				fix_position[32 - m - 1] = false;

		}

		if (invalidInput == true)
		{
			invalidInput = false;
			continue;
		}

		//fprintf(stderr,  "\n ----------------------->kmer_to_search: %s", kmer_to_search);
		//fprintf(stderr,"%d \n", i);
		//Cod and r.c.
		pattern_to_search = create_pattern_to_search(kmer_to_search, 32);

		//Store all information regard the search in seq+
		seq = init_bwt_seq(pattern_to_search, fix_position, 32, max_offset, offset, dp);

		//Core method that solves the approximate pattern matching problem
		get_approximate_match(seq, out, &m_hit_found_for, &m_hit_found_rev);

		free(pattern_to_search);
		free(seq);
	}

	for (int h = out->n_for; h < m_hit_found_for; h++)
		free(out->res_for[h]);

	for (int h = out->n_rev; h < m_hit_found_rev; h++)
		free(out->res_rev[h]);

	free(kmer_to_search);
	free(fix_position);

	//fprintf(stderr,  "\n -------------_>:RETURN OUT");
	return out;
}

/**
 *
 */
void lib_aln_sr_destroy(output* result)
{
	//fprintf(stderr, "result->n_for: %"PRIu16 "-", result->n_for);

	for (int i = 0; i < result->n_for; i++)
	{
		free(result->res_for[i]->positions_to_ref);
		if (result->res_for[i]->altSnpPosInsideKmer != NULL)
			free(result->res_for[i]->altSnpPosInsideKmer);

		if (result->res_for[i]->n_mismatches > 0)
			free(result->res_for[i]->different_positions);

		free(result->res_for[i]);
	}

	for (int i = 0; i < result->n_rev; i++)
	{
		free(result->res_rev[i]->positions_to_ref);
		if (result->res_rev[i]->altSnpPosInsideKmer != NULL)
			free(result->res_rev[i]->altSnpPosInsideKmer);

		if (result->res_rev[i]->n_mismatches > 0)
			free(result->res_rev[i]->different_positions);

		free(result->res_rev[i]);
	}

	free(result->res_for);
	free(result->res_rev);
	free(result);
}

/**
 *
 */
void lib_aln_hit_reassembly_destroy(hit_reassembly* hr)
{
	free(hr->positions_to_ref);
	free(hr);
}

//To improve and only for  GPLv3 version..
void lib_aln_index(const char* path_genome, const char* prefix, int algo_type)
{
	if (path_genome == 0)
	{
		fprintf(stderr, "Miss pattern to search.\n");
		exit(EXIT_FAILURE);
	}
	else if (path_genome == 0)
	{
		fprintf(stderr, "Miss prefix.\n");
		exit(EXIT_FAILURE);
	}
	else if (algo_type < 0 || algo_type > 4)
	{
		fprintf(stderr, "Miss type of algorithm to apply.\n");
		exit(EXIT_FAILURE);
	}

	fmd_idx_build(path_genome, prefix, algo_type, 10000000); //10000000 is block_size
}

/**
 *
 */
def_param* lib_aln_set_default_parameters(const bwaidx_t* idx_reference, const bwaidx_t* idx_reassembly, const SNPSmartDictionary* snpSD,
											const char quality_threshold, const uint8_t max_mismatches)
{
	def_param* dp = (def_param*) malloc(sizeof(def_param));
	dp->idx_reference = idx_reference;
	dp->idx_reassembly = idx_reassembly;
	dp->snpSD = snpSD;
	dp->max_diff = max_mismatches;
	dp->quality_threshold = quality_threshold;
	return dp;
}

