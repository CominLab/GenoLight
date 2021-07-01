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
#include <inttypes.h>
#include <stdbool.h>
#include "bounded_backtracking_seach.h"
#include "sa.h"
#include "occ.h"
#include "bntseq.h"

#define POS_AMBIGUOUS ((uint32_t)(-1))

//Only for debug: 0 no output,3 output
int verbose_bound_backtracking_search = 0;

static int cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width);

static stack_t * init_stack(int nmismatch);
static void reset_stack(stack_t *stack);
static void destroy_stack(stack_t *stack);

static inline void pop(stack_t *stack, entry_t *e);

static inline void push(stack_t *stack, int i, uint64_t k, uint64_t l, int n_mm, int is_diff);
static inline void shadow(int x, int len, uint64_t max, int last_diff_pos, bwt_width_t *w);

static int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, uint64_t *k0, uint64_t *l0);
static uint64_t* get_pos_from_sa_interval(const bwaidx_t* idx, const uint64_t k, const uint64_t l, uint32_t len, uint64_t* numer_forward,
											uint64_t* numer_rc);

//static void get_string_to_pos(const bwaidx_t *idx, input_query* info_seq, search_result* sr, bool is_rev_comp);

static inline void add_entry_to_result(const bwaidx_t* idx, search_result** rs, uint8_t n_mm, uint8_t start, uint64_t n, uint64_t* pos,
								input_query* info_seq, int* aln, uint8_t offset, bool is_rev_comp);

static inline void add_entry_to_result_snp(search_result** rs, uint64_t hit_code, snp_kmer* sk, uint32_t* pos_mismatch, uint8_t n_mm, uint8_t offset,
									int* aln);

static int comp(const void * elem1, const void * elem2);

uint64_t encode_hit(const ubyte_t* seq, int len);

//Only for debug
static void print_pattern_to_search(const ubyte_t * seq, int len);

static void print_character(int i);

static uint64_t rev_compl(const uint64_t orig);

static uint32_t* positions_mis_matches(input_query* info_seq, ubyte_t* hit_code, bool is_rev_comp);

void decode_kmer_vargeno(uint64_t e, uint8_t k_length);

snp_kmer* GetPosFromDictionarySNP(uint8_t isAmbiguos, uint32_t pos, uint8_t altSnpPosInsideKmer, snp_aux_table* snp_aux_tab);

//Modified version of bwt_match_gap in bwtgap.c file
hit_reassembly* get_exact_match_to_compleate_dictionary(bwaidx_t *idx, ubyte_t* seq, int len)
{

	//Max number(s) of mismatch(es) between a hit and the reference
	const int max_diff = 0;

	//Warning: Hit means an admissible match (with respect to the input parameters) between the pattern and the reference

	if (verbose_bound_backtracking_search > 2)
		fprintf(stderr, "Start get_approximate_match\n");

	/*
	 * Number of permissible match found.
	 * Warning: if the current hit found and its r.v. are present inside the reference, it's increased by two
	 */

	int n_hit_found = 0;

	hit_reassembly* res = (hit_reassembly*) malloc(sizeof(hit_reassembly));

	/*
	 * Heap-like data structure to keep partial hits. It is prioritized on the number of mismatch
	 * inside the partial hits(also called entry): less mismatch have a entry first is extract
	 */
	stack_t *stack = init_stack(max_diff);

	//see cal_width implementation
	bwt_width_t* width = (bwt_width_t*) malloc((len + 1) * sizeof(bwt_width_t));
	cal_width(idx->bwt, len, seq, width);

	//Get bwt
	bwt_t* bwt = idx->bwt;

	//Print pattern to search
	if (verbose_bound_backtracking_search > 2)
	{
		fprintf(stderr, "\nThe search pattern is the r.c. of input pattern: ");
		print_pattern_to_search(seq, len);
	}
	//Reset stack
	reset_stack(stack);

	push(stack, len, 0, bwt->seq_len, 0, 0);

	//With priority on the number of mismatches, partial hits calculated previously are extracted
	while (stack->n_entries)
	{
		//Get the best entry among the previously calculated partial hits
		entry_t e;
		pop(stack, &e);

		//Defines whether the current entry 'e' is a hit
		bool hit_found = false;

		/*
		 * The number of mismatch of the current entry regard the substring [i,n-1] of the
		 * input pattern.
		 */
		int i = e.info;

		//(k,l) is the SA region of [i,n-1]
		uint64_t k = e.k;
		uint64_t l = e.l;

		/*
		 * For the current entry 'e' previously calculated, n_mm_yet_allow is the max number of
		 * mismatch that 'e' can still contain
		 */
		int n_mm_yet_allow = 0 - e.n_mm; // 0 = max_diff

		if (verbose_bound_backtracking_search > 2)
		{
			fprintf(stderr, "\n*New entry is extracted, this has: ");
			fprintf(stderr, "\ni: %d\n", i);
			fprintf(stderr, "k: %ld\n", k);
			fprintf(stderr, "l: %ld\n", l);
			fprintf(stderr, "n_mm_yet_allow: %d\n\n", n_mm_yet_allow);
		}

		//This control should be useless..
		if (n_mm_yet_allow < 0)
			continue;

		if (i > 0 && n_mm_yet_allow < width[i - 1].bid)
		{
			if (verbose_bound_backtracking_search > 2)
			{
				fprintf(stderr, "°°Used width!°°\n");
				fprintf(stderr, "n_mm_yet_allow: %d\n", n_mm_yet_allow);
				fprintf(stderr, "width[i - 1].bid: %d\n", width[i - 1].bid);
			}
			continue;
		}

		// Check whether a hit is found
		if (i == 0)
		{
			//This means that the length of the current entry is exactly equal to the length of the input patter
			hit_found = true;
		}
		else if (n_mm_yet_allow == 0)
		{
			/*
			 * This means that the current partial hit have the max possible number of mismatch,
			 * so the only possible hit is composed by 'e' and the miss base must be equal to the
			 * base of the input patter
			 */
			if (verbose_bound_backtracking_search > 2)
				fprintf(stderr, "bwt_match_exact_alt\n");

			if (bwt_match_exact_alt(bwt, i, seq, &k, &l))
			{
				if (verbose_bound_backtracking_search > 2)
				{
					fprintf(stderr, "hit found\n");
					fprintf(stderr, "k: %ld\n", k);
					fprintf(stderr, "l: %ld\n", l);
				}
				hit_found = true;
			}
			else
			{
				if (verbose_bound_backtracking_search > 2)
					fprintf(stderr, "Hit not found from bwt_match_exact_alt\n");

				continue; // no hit, skip
			}
		}
		if (hit_found)
		{
			//This method is used to reduce the search space
			shadow(l - k + 1, len, bwt->seq_len, e.last_diff_pos, width);

			//We might have to add both the direct and the reverse complement of the current entry
			if (n_hit_found == 1)
			{
				fprintf(stderr, "More k-mer found inside re-assembly\n");
				exit(EXIT_FAILURE);
			}

			/*
			 if (verbose_bound_backtracking_search > 2)
			 {
			 printf("k: %"PRIu64 "\n", k);
			 printf("l: %"PRIu64 "\n", l);
			 }
			 */

			//Define if pos is the start position of the input patter inside the reference genome or the r.c of the input pattern
			bool strand;

			//Max number of occurrences of current hit found inside the reference
			int max_occ = l - k + 1;

			//fprintf(stderr, "max_occ: %d\n", max_occ);
			uint64_t t;
			res->positions_to_ref = (uint64_t*) malloc(sizeof(uint64_t));

			for (t = k; t <= l; ++t)
			{
				uint64_t pos = bwa_sa2pos(idx->bns, idx->bwt, t, len, &strand);

				if (pos != ULLONG_MAX)
				{
					if (strand == 0)
					{
						res->positions_to_ref[0] = pos;
						res->isRev = false;
					}
					else if (strand == 1)
					{
						res->positions_to_ref[0] = pos;
						res->isRev = true;
					}
					break;
				}
			}

			continue;
		}

		//Increase(not decrease..) the size of the current partial hit
		--i;

		/*
		 * Give k - 1 and l of the current entry 'e', compute
		 * O(x,k-1) and O(x,l), where x can be A,C,G,T.
		 *
		 * For more detail please see section (2.3) and (2.4) of original paper:
		 * https://academic.oup.com/bioinformatics/article/25/14/1754/225615
		 *
		 */
		uint64_t cnt_k[4], cnt_l[4];
		bwt_2occ4(bwt, k - 1, l, cnt_k, cnt_l);

		//Try to extend the current partial hit whit every possible base
		for (int j = 1; j <= 4; ++j)
		{
			int c = (seq[i] + j) & 3;

			if (verbose_bound_backtracking_search > 2)
			{
				fprintf(stderr, ">Try to extend the current partial hit whit: ");
				print_character(c);
			}
			int is_mm = (j != 4 || seq[i] > 3);
			k = bwt->L2[c] + cnt_k[c] + 1;
			l = bwt->L2[c] + cnt_l[c];

			if (k <= l)
			{
				if (verbose_bound_backtracking_search > 2)
				{
					if (is_mm)
						fprintf(stderr, "Found(whit mismatch)-->push inside stack this partial hit whit this values: \n");
					else
						fprintf(stderr, "Found(no mismatch)-->push inside stack this partial hit whit this values: \n");

					fprintf(stderr, "i %d\n", i);
					fprintf(stderr, "k %ld:\n", k);
					fprintf(stderr, "l %ld:\n", l);
					fprintf(stderr, "num. mismatch inside this partial hit is: %d\n", e.n_mm + is_mm);
				}
				//The partial hit is add in the reference
				push(stack, i, k, l, e.n_mm + is_mm, is_mm);
			}
			else
			{
				if (verbose_bound_backtracking_search > 2)
					fprintf(stderr, "Not found-->discard\n");
			}
		}

	}

	if (verbose_bound_backtracking_search > 2)
	{
		fprintf(stderr, "End exact match!\n");
	}

//free memory
	destroy_stack(stack);
	free(width);

	return res;
}

//Modified version of bwt_match_gap in bwtgap.c file
void get_approximate_match(input_query* info_seq, output* out, uint16_t* m_hit_found_for, uint16_t* m_hit_found_rev)
{

//Warning: Hit means an admissible match (with respect to the input parameters) between the pattern and the reference

	//if (verbose_bound_backtracking_search > 2)
	//fprintf(stderr, "--------------------->Start get_approximate_match\n");

	/*
	 * Number of permissible match found.
	 * Warning: if the current hit found and its r.v. are present inside the reference, it's increased by two
	 */
	int n_hit_found_for = out->n_for;
	int n_hit_found_rev = out->n_rev;

	search_result** res_for = out->res_for;
	search_result** res_rev = out->res_rev;

	bool* fix_pos = info_seq->fix_pos;

	bwaidx_t *idx = NULL;
	bool need_dictionary = false;
	struct Entry_dictionary_snp* dictionary = info_seq->dp->snpSD->snp_dict;

	for (int i = 0; i < 2; i++)
	{
		if (i == 0)
		{
			idx = info_seq->dp->idx_reference;
			need_dictionary = false;
			//fprintf(stderr, "\n------>NO DIZIONARIO******* ");
		}
		else
		{
			idx = info_seq->dp->idx_reassembly;
			need_dictionary = true;
			//fprintf(stderr, "\n----->INIZIO A CERCARE NEL DIZIONARIO***** ");
		}

		/*
		 * Heap-like data structure to keep partial hits. It is prioritized on the number of mismatch
		 * inside the partial hits(also called entry): less mismatch have a entry first is extract
		 */
		stack_t *stack = init_stack(info_seq->dp->max_diff);

		//see cal_width implementation
		bwt_width_t* width = (bwt_width_t*) malloc((info_seq->len + 1) * sizeof(bwt_width_t));

		//Get bwt
		bwt_t* bwt = idx->bwt;

		cal_width(bwt, info_seq->len, info_seq->seq, width);

		//Max number(s) of mismatch(es) between a hit and the reference
		const int max_diff = info_seq->dp->max_diff;

		//Print pattern to search
		if (verbose_bound_backtracking_search > 2)
		{
			fprintf(stderr, "\nThe search pattern is the r.c. of input pattern: ");
			print_pattern_to_search(info_seq->seq, info_seq->len);
		}
		//Reset stack
		reset_stack(stack);

		push(stack, info_seq->len, 0, bwt->seq_len, 0, 0);

		//With priority on the number of mismatches, partial hits calculated previously are extracted
		while (stack->n_entries)
		{
			//Get the best entry among the previously calculated partial hits
			entry_t e;
			pop(stack, &e);

			//Defines whether the current entry 'e' is a hit
			bool hit_found = false;

			/*
			 * The number of mismatch of the current entry regard the substring [i,n-1] of the
			 * input pattern.
			 */
			int i = e.info;

			//(k,l) is the SA region of [i,n-1]
			uint64_t k = e.k;
			uint64_t l = e.l;

			/*
			 * For the current entry 'e' previously calculated, n_mm_yet_allow is the max number of
			 * mismatch that 'e' can still contain
			 */
			int n_mm_yet_allow = max_diff - e.n_mm;

			if (verbose_bound_backtracking_search > 2)
			{
				fprintf(stderr, "\n*New entry is extracted, this has: ");
				fprintf(stderr, "\ni: %d\n", i);
				fprintf(stderr, "k: %ld\n", k);
				fprintf(stderr, "l: %ld\n", l);
				fprintf(stderr, "n_mm_yet_allow: %d\n\n", n_mm_yet_allow);
			}

			//This control should be useless..
			if (n_mm_yet_allow < 0)
				continue;

			if (i > 0 && n_mm_yet_allow < width[i - 1].bid)
			{
				if (verbose_bound_backtracking_search > 2)
				{
					fprintf(stderr, "°°Used width!°°\n");
					fprintf(stderr, "n_mm_yet_allow: %d\n", n_mm_yet_allow);
					fprintf(stderr, "width[i - 1].bid: %d\n", width[i - 1].bid);
				}
				continue;
			}

			// Check whether a hit is found
			if (i == 0)
			{
				//This means that the length of the current entry is exactly equal to the length of the input patter
				hit_found = true;
			}
			else if (n_mm_yet_allow == 0)
			{
				/*
				 * This means that the current partial hit have the max possible number of mismatch,
				 * so the only possible hit is composed by 'e' and the miss base must be equal to the
				 * base of the input patter
				 */
				if (verbose_bound_backtracking_search > 2)
					fprintf(stderr, "bwt_match_exact_alt\n");

				if (bwt_match_exact_alt(bwt, i, info_seq->seq, &k, &l))
				{
					if (verbose_bound_backtracking_search > 2)
					{
						fprintf(stderr, "hit found\n");
						fprintf(stderr, "k: %ld\n", k);
						fprintf(stderr, "l: %ld\n", l);
					}
					hit_found = true;
				}
				else
				{
					if (verbose_bound_backtracking_search > 2)
						fprintf(stderr, "Hit not found from bwt_match_exact_alt\n");

					continue; // no hit, skip
				}
			}
			if (hit_found)
			{
				//This method is used to reduce the search space
				shadow(l - k + 1, info_seq->len, bwt->seq_len, e.last_diff_pos, width);
				//printf("\nm_hit_found_for: %"PRIu16 "\n", (*m_hit_found_for));

				//We might have to add both the direct and the reverse complement of the current entry
				if (n_hit_found_for + 1 > (*m_hit_found_for))
				{
					//fprintf(stderr, "\n°°°°°°°°°°°°°°Start realloc for \n");

					(*m_hit_found_for) <<= 1;
					res_for = (search_result**) realloc(res_for, (*m_hit_found_for) * sizeof(search_result*));

					for (int i = n_hit_found_for; i < (*m_hit_found_for); i++)
						*(res_for + i) = (search_result *) malloc(sizeof(search_result));

					//fprintf(stderr, "°°°°°°°°°°°°°Realloc successfully completed\n");
				}

				if (n_hit_found_rev + 1 > (*m_hit_found_rev))
				{
					//fprintf(stderr, "\n°°°°°°°°°°°°°°Start realloc rev \n");

					(*m_hit_found_rev) <<= 1;
					res_rev = (search_result**) realloc(res_rev, (*m_hit_found_rev) * sizeof(search_result*));

					for (int i = n_hit_found_rev; i < (*m_hit_found_rev); i++)
						*(res_rev + i) = (search_result *) malloc(sizeof(search_result));

					//fprintf(stderr, "\n°°°°°°°°°°°°°°°°°Realloc successfully completed\n");
				}

				/*
				 if (verbose_bound_backtracking_search > 2)
				 {
				 printf("k: %"PRIu64 "\n", k);
				 printf("l: %"PRIu64 "\n", l);
				 }
				 */
				if (need_dictionary)
				{
					//fprintf(stderr, "NEED DICTIONARY!\n");

					//For each index j that belongs to the suffix interval [k,l], get_pos_from_sa_interval
					//return all SA(j) valid

					//Variable that stores if info_seq->seq,  have the same orientation of hit_code
					bool is_kmerToSearch_rc_respect_reassembly;

					uint64_t position_to_reassembly = ULLONG_MAX;
					bool strand;
					//Max number of occurrences of current hit found inside the reference
					int max_occ = l - k + 1;

					//fprintf(stderr, "max_occ: %d\n", max_occ);
					uint64_t t;

					for (t = k; t <= l; ++t)
					{
						position_to_reassembly = bwa_sa2pos(idx->bns, idx->bwt, t, 32, &strand);

						if (position_to_reassembly != ULLONG_MAX)
						{
							if (strand == 0)
								is_kmerToSearch_rc_respect_reassembly = false;
							else
								is_kmerToSearch_rc_respect_reassembly = true;
							break;
						}
					}

					//If all positions corresponding to suffix interval [k, l] are not valid continue
					if (position_to_reassembly == ULLONG_MAX)
						continue;

					//fprintf(stderr, "pos to reassembly: %"PRIu64 "\n", position_to_reassembly);

					//Nel caso del dizionario degli snp ha un solo elemento
					struct Entry_dictionary_snp* entryDictsnp = &dictionary[position_to_reassembly];

					/*
					value_snp* entryDictsnp = principalDictionary_get(info_seq->dp->snpSD->snp_dict, position_to_reassembly);
					if (entryDictsnp == NULL)
											continue;
					*/
					if ((entryDictsnp->pos_snp_genome == 0) && (entryDictsnp->isDirectAmb == 0) && (entryDictsnp->isRevAmb == 0))
					{
						//I found a spourius k-mer between two contigs
						continue;
					}

					//fprintf(stderr, "Bit rev: %" PRIu8 "\n", entryDictsnp.is_reverse);
					//fprintf(stderr, "Is_compleate: %" PRIu8 "\n", entryDictsnp.is_compleate);

					if (entryDictsnp->pos_snp_genome == POS_AMBIGUOUS)
						continue;

					int64_t len = 0;

					//Get string that start in position_to_reassembly[0] inside reassembly
					ubyte_t* hit_code = bns_get_seq(idx->bns->l_pac, idx->pac, position_to_reassembly, position_to_reassembly + info_seq->len, &len);

					if (info_seq->len != len)			//This must false..
					{
						printf("Impossible to extract the patter from the reference\n");
						exit(EXIT_FAILURE);
					}

					//Get snp_auxiliary_table
					snp_aux_table* snp_aux_tab = info_seq->dp->snpSD->snp_auxiliary_table;

					uint32_t* pos_mis_match = NULL;

					if (e.n_mm > 0)
						pos_mis_match = positions_mis_matches(info_seq, hit_code, is_kmerToSearch_rc_respect_reassembly);

					snp_kmer* sk;

					//fprintf(stderr, "Bit rev: %" PRIu8 "\n", entryDictsnp.is_reverse);
					//fprintf(stderr, "Is_compleate: %" PRIu8 "\n", entryDictsnp.is_compleate);

					/*
					 * If the i-th k-mer extracted from the current reads, stored inside "info_seq->seq", is present whit the same orientation
					 * inside reference genome(modified with snp), means that:
					 *
					 * -Is found no r.c. inside re-assembly(has_kmer_same_orientation_reassembly = false), and inside re-assembly the k-mer is the same that
					 * 	store inside reference genome(entryDictsnp.is_reverse = false)
					 *
					 * -Is found r.c. inside reassembly(has_kmer_same_orientation_reassembly = true), but the k-mer inside reassembly is the r.c. of the k-mer
					 *  present inside reference genome(entryDictsnp.is_reverse = true)
					 *
					 * So, if (entryDictsnp.is_reverse == has_kmer_same_orientation_reassembly) means that the k-mer extracted from the reads is present
					 * exactly in the reference genome(modified with snp)
					 */

					if ((entryDictsnp->is_compleate == 0) && (entryDictsnp->is_reverse == is_kmerToSearch_rc_respect_reassembly))
					{
						//In case it's direct
						if (is_kmerToSearch_rc_respect_reassembly == false)
						{
							//fprintf(stderr, "614\n");
							//pos_mis_math è già giusta

							sk = GetPosFromDictionarySNP(entryDictsnp->isDirectAmb, entryDictsnp->pos_snp_genome, entryDictsnp->altSnpPosInsideKmer,
															snp_aux_tab);

							add_entry_to_result_snp(res_for, encode_hit(hit_code, 32), sk, pos_mis_match, e.n_mm, info_seq->offset, &n_hit_found_for);

						}
						else				//This means has_kmer_same_orientation_reassembly == entryDictsnp.is_reverse == true
						{

							//fprintf(stderr, "620\n");

							sk = GetPosFromDictionarySNP(entryDictsnp->isRevAmb, entryDictsnp->pos_snp_genome, entryDictsnp->altSnpPosInsideKmer,
															snp_aux_tab);
							// In this case I have to do the r.c. of the k-mer found in the reassembly(hit_code) as
							// it is the r.c of the searched one (info_seq->seq) and of the one present in the reference genome (entryDictsnp.is_reverse = true)

							add_entry_to_result_snp(res_for, rev_compl(encode_hit(hit_code, 32)), sk, pos_mis_match, e.n_mm, info_seq->offset,
													&n_hit_found_for);
						}
					}

					/*
					 * If the current k-mer isn't complete and presents different orientation  of that stored in the smart dictionary in memory means that the current
					 *  k-mer isn't present inside the reference genome but is present his r.c
					 */

					/*
					 * If the i-th k-mer extracted from the current reads, stored inside "info_seq->seq", is present r.c. inside reference genome(modified with snp),
					 * means that:
					 *
					 * -Is found no r.c. inside re-assembly(has_kmer_same_orientation_reassembly = false), and inside re-assembly the k-mer is the r.c.
					 *  respect k-mer inside reference genome(entryDictsnp.is_reverse = true)
					 *
					 * -Is found r.c. inside reassembly(has_kmer_same_orientation_reassembly = true), but the k-mer inside reassembly is the same of the k-mer present
					 * inside reference genome(entryDictsnp.is_reverse = false)
					 *
					 * So, if (entryDictsnp.is_reverse != has_kmer_same_orientation_reassembly) means that the k-mer extracted from the reads is present
					 * r.c. in the reference genome(modified with snp)
					 */

					else if ((entryDictsnp->is_compleate == 0) && (entryDictsnp->is_reverse != is_kmerToSearch_rc_respect_reassembly))
					{
						if (is_kmerToSearch_rc_respect_reassembly == false)
						{
							//fprintf(stderr, "659\n");

							//if (entryDictsnp.num_occurence_direct == 0)
							//	fprintf(stderr, "prima di entrare in GetPosFrom stampo il valore: %"PRIu32 "\n", entryDictsnp.pos_snp_genome);

							if (e.n_mm > 0)
								pos_mis_match[0] = info_seq->len - pos_mis_match[0] + 1;

							sk = GetPosFromDictionarySNP(entryDictsnp->isDirectAmb, entryDictsnp->pos_snp_genome, entryDictsnp->altSnpPosInsideKmer,
															snp_aux_tab);

							//if (entryDictsnp.num_occurence_direct == 0)
							//      fprintf(stderr, "real_pos_ref[0]: %"PRIu32 "\n", real_pos_ref[0];

							add_entry_to_result_snp(res_rev, rev_compl(encode_hit(hit_code, 32)), sk, pos_mis_match, e.n_mm,
													info_seq->max_offset - info_seq->offset, &n_hit_found_rev);
						}
						else
						{
							//fprintf(stderr, "670\n");
							if (e.n_mm > 0)
								pos_mis_match[0] = info_seq->len - pos_mis_match[0] + 1;

							sk = GetPosFromDictionarySNP(entryDictsnp->isRevAmb, entryDictsnp->pos_snp_genome, entryDictsnp->altSnpPosInsideKmer,
															snp_aux_tab);
							add_entry_to_result_snp(res_rev, encode_hit(hit_code, 32), sk, pos_mis_match, e.n_mm,
													info_seq->max_offset - info_seq->offset, &n_hit_found_rev);
						}
					}

					/*
					 * Remember:
					 * If current k-mer is "complete", means that it's present both directly and r.c. in reference genome.
					 *
					 * Since all the reads they could be processed first directly and then r.c., store information of both the current k-mer
					 * and his r.c inside,respectively, res_for and res_rev
					 */

					else if (entryDictsnp->is_compleate == 1)
					{
						value* element = total_kmer_snp_struct_get(info_seq->dp->snpSD->reverse_total_kmer_snp, position_to_reassembly);

						//fprintf(stderr, "700\n");

						if (element == NULL)
						{
							fprintf(stderr, "Error: reference return by unordered map(for k-mer compleate) is NULL\n");
							exit(EXIT_FAILURE);
						}

						uint32_t pos_unord_map = element->pos;

						uint8_t altSnpPosRC = element->snpInfo;

						uint32_t* pos_mis_matchRev = NULL;

						//fprintf(stderr, "pos inside unordered map: %" PRIu32 "\n", pos_unord_map);
						//fprintf(stderr, "position_to_reassembly: %" PRIu32 "\n",position_to_reassembly );

						// Always in the "main" dictionary I have the info of the k-mer that appears in the reassembly while in the unordered map its r.c.

						if (is_kmerToSearch_rc_respect_reassembly == false)	//if the k-mer search has the same orientation as the k-mer stored in the reassembly
						{

							//fprintf(stderr, "712\n");

							//fprintf(stderr, "714 INDEX: %"PRIu32 "\n", entryDictsnp.pos_snp_genome);
							//Memorizzo
							sk = GetPosFromDictionarySNP(entryDictsnp->isDirectAmb, entryDictsnp->pos_snp_genome, entryDictsnp->altSnpPosInsideKmer,
															snp_aux_tab);
							add_entry_to_result_snp(res_for, encode_hit(hit_code, 32), sk, pos_mis_match, e.n_mm, info_seq->offset, &n_hit_found_for);

							if (pos_unord_map != POS_AMBIGUOUS)
							{
								if (e.n_mm > 0)
								{
									pos_mis_matchRev = (uint32_t*) malloc(sizeof(uint32_t));
									pos_mis_matchRev[0] = info_seq->len - pos_mis_match[0] + 1;
									//fprintf(stderr, "pos mis_match: %" PRIu32 "\n", pos_mis_match[0]);
									//fprintf(stderr, "pos_mis_matchRev: %" PRIu32 "\n", pos_mis_matchRev[0]);
								}

								snp_kmer* skRev = GetPosFromDictionarySNP(entryDictsnp->isRevAmb, pos_unord_map, altSnpPosRC, snp_aux_tab);

								add_entry_to_result_snp(res_rev, rev_compl(encode_hit(hit_code, 32)), skRev, pos_mis_matchRev, e.n_mm,
														info_seq->max_offset - info_seq->offset, &n_hit_found_rev);
							}
						}
						//if the k-mer I'm processing has a different orientation than the k-mer stored in the dictionary
						else	//has_kmer_same_orientation_reassembly == true
						{
							//fprintf(stderr, "730\n");

							if (pos_unord_map != POS_AMBIGUOUS)
							{

								//fprintf(stderr, "749\n");
								//fprintf(stderr, "position_to_reassembly: %"PRIu64 "\n",position_to_reassembly);
								//fprintf(stderr, "pos_unord_map: %"PRIu32 "\n",pos_unord_map);

								sk = GetPosFromDictionarySNP(entryDictsnp->isRevAmb, pos_unord_map, altSnpPosRC, snp_aux_tab);
								add_entry_to_result_snp(res_for, rev_compl(encode_hit(hit_code, 32)), sk, pos_mis_match, e.n_mm, info_seq->offset,
														&n_hit_found_for);
							}
							//fprintf(stderr, "position inside reference genoma:: %"PRIu32 "\n",entryDictsnp.pos_snp_genome );

							if (e.n_mm > 0)
							{
								pos_mis_matchRev = (uint32_t*) malloc(sizeof(uint32_t));
								pos_mis_matchRev[0] = info_seq->len - pos_mis_match[0] + 1;
							}

							snp_kmer* skRev = GetPosFromDictionarySNP(entryDictsnp->isDirectAmb, entryDictsnp->pos_snp_genome,
																		entryDictsnp->altSnpPosInsideKmer, snp_aux_tab);
							//fprintf(stderr, "759\n");
							add_entry_to_result_snp(res_rev, encode_hit(hit_code, 32), skRev, pos_mis_matchRev, e.n_mm,
													info_seq->max_offset - info_seq->offset, &n_hit_found_rev);
						}

					}
					free(hit_code);
				}

				else
				{

					/*
					 * Inside the suffix interval [k,l], some position is relative
					 * of a hit that is forward respect the pattern input, other
					 * is r.c. respect the pattern input
					 */
					uint64_t numer_forward = 0;
					uint64_t numer_revC = 0;

					//For each index j that belongs to the suffix interval [k,l], get_pos_from_sa_interval
					//return all SA(j) valid
					uint64_t* position_to_reference_genome = get_pos_from_sa_interval(idx, k, l, info_seq->len, &numer_forward, &numer_revC);

					//If all positions corresponding to suffix interval [k, l] are not valid continue
					if (((numer_forward + numer_revC) == 0) || (position_to_reference_genome == NULL))
						continue;

					if (numer_forward > 0)
					{
						//fprintf(stderr, "\n 791 \n");
						add_entry_to_result(idx, res_for, e.n_mm, 0, numer_forward, position_to_reference_genome, info_seq, &n_hit_found_for,
											info_seq->offset, false);
					}
					if (numer_revC > 0)
					{
						//fprintf(stderr, "\797\n");
						add_entry_to_result(idx, res_rev, e.n_mm, numer_forward, numer_revC, position_to_reference_genome, info_seq, &n_hit_found_rev,
											info_seq->max_offset - info_seq->offset, true);
					}
					free(position_to_reference_genome);
				}
				continue;
			}

			//Increase(not decrease..) the size of the current partial hit
			--i;

			/*
			 * Give k - 1 and l of the current entry 'e', compute
			 * O(x,k-1) and O(x,l), where x can be A,C,G,T.
			 *
			 * For more detail please see section (2.3) and (2.4) of original paper:
			 * https://academic.oup.com/bioinformatics/article/25/14/1754/225615
			 *
			 */
			uint64_t cnt_k[4], cnt_l[4];
			bwt_2occ4(bwt, k - 1, l, cnt_k, cnt_l);

			if ((fix_pos[i] == false) || (fix_pos == NULL))
			{
				//Try to extend the current partial hit whit every possible base
				for (int j = 1; j <= 4; ++j)
				{
					int c = (info_seq->seq[i] + j) & 3;

					if (verbose_bound_backtracking_search > 2)
					{
						fprintf(stderr, ">Try to extend the current partial hit whit: ");
						print_character(c);
					}
					int is_mm = (j != 4 || info_seq->seq[i] > 3);
					k = bwt->L2[c] + cnt_k[c] + 1;
					l = bwt->L2[c] + cnt_l[c];

					if (k <= l)
					{
						if (verbose_bound_backtracking_search > 2)
						{
							if (is_mm)
								fprintf(stderr, "Found(whit mismatch)-->push inside stack this partial hit whit this values: \n");
							else
								fprintf(stderr, "Found(no mismatch)-->push inside stack this partial hit whit this values: \n");

							fprintf(stderr, "i %d\n", i);
							fprintf(stderr, "k %ld:\n", k);
							fprintf(stderr, "l %ld:\n", l);
							fprintf(stderr, "num. mismatch inside this partial hit is: %d\n", e.n_mm + is_mm);
						}
						//The partial hit is add in the reference
						push(stack, i, k, l, e.n_mm + is_mm, is_mm);
					}
					else
					{
						if (verbose_bound_backtracking_search > 2)
							fprintf(stderr, "Not found-->discard\n");
					}
				}
			}
			else
			{
				//i corresponds to a position in which character can not be changed
				int c = info_seq->seq[i];
				int is_mm = 0;
				k = bwt->L2[c] + cnt_k[c] + 1;
				l = bwt->L2[c] + cnt_l[c];

				if (k <= l)
				{
					push(stack, i, k, l, e.n_mm + is_mm, is_mm);
				}
			}
		}

		/*
		 printf("\nm_hit_found_for: %"PRIu16 "\n", (*m_hit_found_for));
		 printf("m_hit_found_rev: %"PRIu16 "\n", (*m_hit_found_rev));

		 printf("\nn_hit_found_for: %"PRIu16 "\n", n_hit_found_for);
		 printf("n_hit_found_rev: %"PRIu16 "\n", n_hit_found_rev);
		 */

		if (verbose_bound_backtracking_search > 2)
		{
			fprintf(stderr, "\nNumber of hit found: %d\n", n_hit_found_for + n_hit_found_rev);
			fprintf(stderr, "End get_approximate_match\n");
		}

		//free memory
		destroy_stack(stack);
		free(width);
	}
	out->n_for = n_hit_found_for;
	out->n_rev = n_hit_found_rev;

	out->res_for = res_for;
	out->res_rev = res_rev;
	//fprintf(stderr, " ---------------------->854-fine get_approximate_match\n");
}

input_query* init_bwt_seq(uint8_t* pattern_to_search, bool* fix_pos, const size_t pattern_len, const uint8_t max_offset, const uint8_t offset,
							const def_param* dp)
{
	input_query *seq = (input_query*) calloc(1, sizeof(input_query));
	seq->seq = pattern_to_search;
	seq->fix_pos = fix_pos;
	seq->len = pattern_len;
	seq->max_offset = max_offset;
	seq->dp = dp;
	seq->offset = offset;
	return seq;
}

uint8_t* create_pattern_to_search(const char* pattern_input, const size_t pattern_len)
{
//Code pattern input
	uint8_t* pattern_to_search = (uint8_t*) malloc(pattern_len * sizeof(uint8_t));

	uint8_t base_to_int[] =

	{0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0xFF, 0x01, 0xFF, 0xFF,
		0xFF, 0x02, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x04, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x03, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0xFF, 0x01, 0xFF, 0xFF, 0xFF, 0x02, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x04, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0x03, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF};

	for (int j = 0; j < pattern_len; ++j)
	{
		pattern_to_search[j] = base_to_int[(size_t) pattern_input[pattern_len - j - 1]]; //reverse
		pattern_to_search[j] = pattern_to_search[j] > 3 ? 4 : 3 - pattern_to_search[j]; //complement
	}
	return pattern_to_search;
}

/*
 *Warning: this method works on the reverse (not complement) of the patter T given by the user.
 *
 *width[i].bid represents the lower bound of the number of differences in P_r[0,i], where
 *P_r is the reverse of the input pattern. This is useful for making the search space smaller.
 */
//Derived from bwt_cal_width in bwtaln.c file
static int cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width)
{
	if (verbose_bound_backtracking_search > 2)
		fprintf(stderr, "\n->Inside bwt_cal_width\n");

//[k,l] is the S.A interval values
	uint64_t k, l;

//ok e ol represent represent the values of the rank function
	uint64_t ok, ol;

	int i, bid;
	bid = 0;

	k = 0;
	l = bwt->seq_len;

	for (i = 0; i < len; ++i)
	{
		//Before i have complement
		ubyte_t c = 3 - str[i];

		if (verbose_bound_backtracking_search > 2)
			fprintf(stderr, "Current character considered: %c\n", "ACGTN"[c]);

		if (c < 4)
		{
			//This method calculates the rank function(usually indicated with O) of c
			bwt_2occ(bwt, k - 1, l, c, &ok, &ol);

			//L2[c] is the Count function (usually indicated with C)
			k = bwt->L2[c] + ok + 1;
			l = bwt->L2[c] + ol;

			if (verbose_bound_backtracking_search > 2)
			{
				fprintf(stderr, "k: %ld\n", k);
				fprintf(stderr, "l: %ld\n", l);
			}
		}
		if (k > l || c > 3)
		{
			// then restart
			if (verbose_bound_backtracking_search > 2)
				fprintf(stderr, "Restart\n");

			k = 0;
			l = bwt->seq_len;
			++bid;
		}
		width[i].w = l - k + 1;
		width[i].bid = bid;

		if (verbose_bound_backtracking_search > 2)
		{
			fprintf(stderr, "--->width[%d].w = %ld\n", i, l - k + 1);
			fprintf(stderr, "--->width[%d].bid = %d\n", i, bid);
		}
	}

	width[len].w = 0;
	width[len].bid = ++bid;
	return bid;
}

//It's the same of gap_reset_stack in bwtgap.c file
static void reset_stack(stack_t *stack)
{
	int i;
	for (i = 0; i != stack->n_sub_stacks; ++i)
		stack->stacks[i].n_entries_substack = 0;
	stack->best = stack->n_sub_stacks;
	stack->n_entries = 0;
}

//Derived from init_stack2 in bwtgap.c file
static stack_t * init_stack(int nmismatch)
{
	stack_t *stack = (stack_t*) calloc(1, sizeof(stack_t));

//Each partial hit, based on the number of mismatch it contains, will be "clustered together",
//in substack(so we can have at most nmismatch + 1 substack)
	stack->n_sub_stacks = nmismatch + 1;
	stack->stacks = (substack_t*) calloc(stack->n_sub_stacks, sizeof(stack_t));
	return stack;
}

//Derived from gap_push in bwtgap.c file
static inline void push(stack_t *stack, int i, uint64_t k, uint64_t l, int n_mm, int is_diff)
{
//Get pointer to substack relative to partial hit whith n_mm mismatch
	substack_t *q = stack->stacks + n_mm;

	if (q->n_entries_substack == q->m_entries)
	{
		//fprintf(stderr, "q->n_entries_substack: %d\n", q->m_entries);
		//fprintf(stderr, "q->m_entries: %d\n", q->m_entries);

		//fprintf(stderr, "Prima: q->m_entries: %d\n", q->m_entries);
		q->m_entries = q->m_entries ? q->m_entries << 1 : 4;
		q->start_substack = (entry_t*) realloc(q->start_substack, sizeof(entry_t) * q->m_entries);
		//fprintf(stderr, "Dopo: q->m_entries: %d\n", q->m_entries);
	}

//Add entry at the end of the substack
	entry_t *p = q->start_substack + q->n_entries_substack;
	p->info = (uint32_t) i;
	p->k = k;
	p->l = l;
	p->n_mm = n_mm;
	p->last_diff_pos = is_diff ? i : 0;

//Increase the total number of entries whit n_mm mismatches
	++(q->n_entries_substack);

//Increase the total number of entries in general
	++(stack->n_entries);

	if (stack->best > n_mm)
	{
		if (verbose_bound_backtracking_search > 2)
			fprintf(stderr, "Update stack-> best\n");
		stack->best = n_mm;
	}
}

//Derived from gap_pop in bwtgap.c file
static inline void pop(stack_t *stack, entry_t *e)
{
	/*
	 * Substack containing partial hits whit number of mismatch is minimum,
	 * all partial hit that NOT belong to this memory region have more mismatch.
	 */
	substack_t *q = stack->stacks + stack->best;

//Of this substack, we take the last one memorized
	*e = q->start_substack[q->n_entries_substack - 1];

//Decrease the total number of entries whit n_mm mismatches
	--(q->n_entries_substack);
//Decrease the total number of entries in general
	--(stack->n_entries);

	/*
	 * If the number of partial hit whit n_mm mismatches now is zero but
	 * the number of entry is not zero, update best for next pop operation
	 */
	if (q->n_entries_substack == 0 && stack->n_entries)
	{
		// reset best
		int i;
		for (i = stack->best + 1; i < stack->n_sub_stacks; ++i)
			if (stack->stacks[i].n_entries_substack != 0)
				break;
		stack->best = i;
	}
	else if (stack->n_entries == 0)
	{
		//Necessary for the first iteration
		stack->best = stack->n_sub_stacks;
	}
}

//It's the same of original function gap_shadow present in bwtgap.c
static inline void shadow(int x, int len, uint64_t max, int last_diff_pos, bwt_width_t *w)
{
	if (verbose_bound_backtracking_search > 2)
	{
		fprintf(stderr, "->Inside shadow\n");
		fprintf(stderr, "last_diff_pos: %d\n", last_diff_pos);
	}
	int i, j;
	for (i = j = 0; i < last_diff_pos; ++i)
	{
		if (verbose_bound_backtracking_search > 2)
		{
			fprintf(stderr, "Current values:\n");
			fprintf(stderr, "w[%d].w: %ld\n", i, w[i].w);
			fprintf(stderr, "w[%d].bid: %d\n", i, w[i].bid);
			fprintf(stderr, "x = l-k+1: %d\n", x);
		}
		if (w[i].w > x)
		{
			w[i].w -= x;
			if (verbose_bound_backtracking_search > 2)
				fprintf(stderr, "Update: w[%d] = %ld: ", i, w[i].w);
		}
		else if (w[i].w == x)
		{
			w[i].bid = 1;
			w[i].w = max - (++j);
			if (verbose_bound_backtracking_search > 2)
			{
				fprintf(stderr, "\nUpdate bid (1)\n");
				fprintf(stderr, "Update w[%d].w: %ld\n", i, w[i].w);
			}
		} // else should not happen
	}

	if (verbose_bound_backtracking_search > 2)
		fprintf(stderr, "->Exit gap_shadow\n");
}

//It's the same of original, present in bwt.c
static int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, uint64_t *k0, uint64_t *l0)
{
	int i;
	uint64_t k, l, ok, ol;
	k = *k0;
	l = *l0;
	for (i = len - 1; i >= 0; --i)
	{
		ubyte_t c = str[i];
		if (c > 3)
			return 0; // there is an N here. no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l)
			return 0; // no match
	}
	*k0 = k;
	*l0 = l;
	return l - k + 1;
}

static uint64_t* get_pos_from_sa_interval(const bwaidx_t* idx, uint64_t k, uint64_t l, uint32_t len, uint64_t* numer_forward, uint64_t* numer_rc)
{

//Current index of the SA, k<=t<=l
	uint64_t t;

//Is SA[t]
	uint64_t pos;

//Define if pos is the start position of the input patter inside the reference genome or the r.c of the input pattern
	bool strand;

//Max number of occurrences of current hit found inside the reference
	int max_occ = l - k + 1;

	if (max_occ > 20)
		return NULL;

//Same meaning of numer_forward and numer_rc
	uint64_t n_f = 0;
	uint64_t n_rc = 0;

	/*
	 * Store all valid positions. The first n_f positions are relative to the forward of the current hit,
	 * the last n_rc to the r.c.
	 */
	uint64_t* set_pos = (uint64_t*) malloc(max_occ * sizeof(uint64_t));
	/*
	 if (verbose_bound_backtracking_search > 2)
	 fprintf(stderr, "Max_occ: %" PRIu64 "\n", max_occ);
	 */
	for (t = k; t <= l; ++t)
	{
		//Calculate SA[t](base 1)
		pos = bwa_sa2pos(idx->bns, idx->bwt, t, len, &strand);

		/*
		 if (verbose_bound_backtracking_search > 2)
		 {
		 if (pos == ULLONG_MAX)
		 fprintf(stderr, "Position return by bwa_sa2pos is not valid\n");
		 else
		 {
		 fprintf(stderr, "Position return by bwa_sa2pos: %" PRIu64 "\n", pos);
		 fprintf(stderr, "Strand: %" PRIu8 "\n", strand);
		 }
		 }
		 */
		//Check if pos is admissible
		if (pos != ULLONG_MAX)
		{
			if (strand == 0)
			{
				set_pos[n_f] = pos;
				n_f++;
			}
			else if (strand == 1)
			{
				set_pos[max_occ - n_rc - 1] = pos;
				n_rc++;
			}
		}
	}

//Number of valid positions found
	uint64_t realSize = n_f + n_rc;

	/*
	 if (verbose_bound_backtracking_search > 2)
	 {
	 fprintf(stderr, "\nn_rc: %"PRIu64 "\n", n_rc);
	 fprintf(stderr, "n_f: %" PRIu64"\n", n_f);
	 }
	 */
	if ((n_rc == n_f) && (n_f == 0))
	{
		if (verbose_bound_backtracking_search > 2)
			fprintf(stderr, "No hits have been added\n");
		return 0;
	}
//Compact set_pos
	if (n_f + n_rc < max_occ)
	{
		for (int i = 0; i < max_occ - realSize; i++)
			set_pos[n_f + i] = set_pos[max_occ - i - 1];
		set_pos = realloc(set_pos, realSize * sizeof(uint64_t));
	}
	/*
	 if (verbose_bound_backtracking_search > 2)
	 {
	 fprintf(stderr, "The positions store inside set_pos are: \n");
	 for (int i = 0; i < realSize; i++)
	 printf("%" PRIu64 "\n", set_pos[i]);
	 }
	 */
	*(numer_forward) = n_f;
	*(numer_rc) = n_rc;

	return set_pos;
}

//Same of gap_destroy_stack in bwtgap.c
static void destroy_stack(stack_t *stack)
{
	int i;
	for (i = 0; i != stack->n_sub_stacks; ++i)
		free(stack->stacks[i].start_substack);
	free(stack->stacks);
	free(stack);
}

//add_entry_to_result(idx, res_rev, e.n_mm, numer_revC, position_to_ref, info_seq, &n_hit_found_rev);

//In the next version of software create a method for add a entry into rs.
inline void add_entry_to_result(const bwaidx_t* idx, search_result** rs, uint8_t n_mm, uint8_t start, uint64_t n, uint64_t* pos, input_query* info_seq,
							int* aln, uint8_t offset, bool is_rev_comp)
{
	if (verbose_bound_backtracking_search > 2)
		fprintf(stderr, "-> Inside add_entry_to_result\n");

	rs[*(aln)]->positions_to_ref = (uint32_t*) malloc(n * sizeof(uint32_t));

	for (int i = 0; i < n; i++)
		rs[*(aln)]->positions_to_ref[i] = pos[i + start];

	rs[*(aln)]->altSnpPosInsideKmer = NULL;
	rs[*(aln)]->offset = offset;
	rs[*(aln)]->num_occur = n;
	rs[*(aln)]->n_mismatches = n_mm;

	uint32_t* different_positions;

	int64_t len = 0;

//Get string that start in position position_to_ref[0]
	ubyte_t* hit_code = bns_get_seq(idx->bns->l_pac, idx->pac, rs[*(aln)]->positions_to_ref[0], rs[*(aln)]->positions_to_ref[0] + info_seq->len,
									&len);

//printf("\n1293-HO TROVATO(ref):\n");
//print_pattern_to_search(hit_code, 32);

	if (info_seq->len != len)			//This must false..
	{
		printf("Impossible to extract the patter from the reference\n");
		exit(EXIT_FAILURE);
	}
	if (n_mm > 0)
	{
		int num_different_pos_found = 0;
		rs[*(aln)]->different_positions = malloc(n_mm * sizeof(uint32_t));

		//Check where the mismatches are
		for (uint32_t j = 0; j < info_seq->len; ++j)
		{
			if ((is_rev_comp == true) && (3 - hit_code[info_seq->len - j - 1] != 3 - info_seq->seq[info_seq->len - j - 1]))
			{
				//fprintf(stderr, "953\n");
				rs[*(aln)]->different_positions[num_different_pos_found] = j + 1; //base-1
				num_different_pos_found++;
			}
			else if ((is_rev_comp == false) && (hit_code[j] != 3 - info_seq->seq[info_seq->len - j - 1]))
			{
				//fprintf(stderr, "959\n");
				rs[*(aln)]->different_positions[num_different_pos_found] = j + 1; //base-1
				num_different_pos_found++;
			}

			if(num_different_pos_found > 1)
			{
				fprintf(stderr, "Error: Impossible num of mis match > 1\n");
				exit(EXIT_FAILURE);
			}
		}
		//fprintf(stderr, "------------------------------------------------\n");

	}
	else
		rs[*(aln)]->different_positions = NULL;

	rs[*(aln)]->hit = encode_hit(hit_code, 32);
	*(aln) = *(aln) + 1;

	free(hit_code);
}

static inline void add_entry_to_result_snp(search_result** rs, uint64_t hit_code, snp_kmer* sk, uint32_t* position_mismatches, uint8_t n_mm, uint8_t offset,
									int* aln)
{
	if (verbose_bound_backtracking_search > 2)
		fprintf(stderr, "-> Inside add_entry_to_result_SNP\n");

	rs[*(aln)]->hit = hit_code;
	rs[*(aln)]->altSnpPosInsideKmer = sk->posAlternativeSnp;
	rs[*(aln)]->positions_to_ref = sk->positions_to_ref;
	rs[*(aln)]->num_occur = sk->num_occur;
	rs[*(aln)]->different_positions = position_mismatches;

	//fprintf(stderr, "n occurance insert IN RES: %" PRIu64 "\n", n_pos_inside_ref_genome);

	rs[*(aln)]->n_mismatches = n_mm;
	rs[*(aln)]->offset = offset;
	*(aln) = *(aln) + 1;

	free(sk);
}

static uint32_t* positions_mis_matches(input_query* info_seq, ubyte_t* hit_code, bool is_rev_comp)
{
	uint8_t num_different_pos_found = 0;

	uint32_t* different_positions = (uint32_t*) malloc(1 * sizeof(uint32_t)); // 1 = n_mm

//Check where the mismatches are
	for (uint32_t j = 0; j < info_seq->len; ++j)
	{
		if ((is_rev_comp == true) && (hit_code[j] != info_seq->seq[j]))
		{
			different_positions[num_different_pos_found] = info_seq->len - j; //base-1
			num_different_pos_found++;
		}
		else if ((is_rev_comp == false) && (hit_code[j] != 3 - info_seq->seq[info_seq->len - j - 1]))
		{
			different_positions[num_different_pos_found] = j + 1; //base-1
			num_different_pos_found++;
		}
	}
	if (num_different_pos_found > 1)
	{
		fprintf(stderr, "Num mm found %" PRIu8 "\n", num_different_pos_found );
		exit(0);
	}
	return different_positions;
}

static void print_pattern_to_search(const ubyte_t * seq, int len)
{
	for (int j = 0; j < len; j++)
		fprintf(stderr, "%c", "ACGTN"[seq[j]]);

	fprintf(stderr, "\n");
}

static void print_character(int i)
{
	char t;

	switch (i)
	{
		case 0:
			t = 'A';
			break;
		case 1:
			t = 'C';
			break;
		case 2:
			t = 'G';
			break;
		case 3:
			t = 'T';
			break;

	}
	fprintf(stderr, "%c\n", t);
}

static int comp(const void * elem1, const void * elem2)
{
	if (*(uint64_t*) elem1 == *(uint64_t*) elem2)
		return 0;
	return *(uint64_t*) elem1 < *(uint64_t*) elem2 ? -1 : 1;
}

uint64_t encode_hit(const ubyte_t* seq, int len)
{
#define KMER_ADD_BASE(x) (encoded_kmer |= (x))

	char * base;
	uint64_t encoded_kmer = 0UL;
	for (int j = 0; j < len; j++)
	{
		encoded_kmer <<= 2;
		switch ("ACGTN"[seq[j]])
		{
			case 'A':
			case 'a':
				KMER_ADD_BASE(0UL);
				//fprintf(stderr, "	A\n");
				break;
			case 'C':
			case 'c':
				KMER_ADD_BASE(1UL);
				//fprintf(stderr, "C\n");
				break;
			case 'G':
			case 'g':
				KMER_ADD_BASE(2UL);
				//fprintf(stderr, "G\n");
				break;
			case 'T':
			case 't':
				KMER_ADD_BASE(3UL);
				//fprintf(stderr, "T\n");
				break;
				break;
		}
	}
#undef KMER_ADD_BASE
	return encoded_kmer;
}

uint64_t rev_compl(const uint64_t orig)
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

	uint64_t ans;
	uint16_t *c_orig = (uint16_t*) (&orig);
	uint16_t *c_ans = (uint16_t*) (&ans);

	for (unsigned i = 0; i < 4; i++)
	{
		c_ans[i] = table[c_orig[3 - i]];
	}

	return ans;
}

void decode_kmer_vargeno(uint64_t e, uint8_t k_length)
{
	char string[32];
	int i;
	for (i = 0; i < k_length; ++i)
	{
		switch (e & 0x0000000000000003)
		{
			case 0:
				string[i] = 'A';
				break;

			case 1:
				string[i] = 'C';
				break;

			case 2:
				string[i] = 'G';
				break;

			case 3:
				string[i] = 'T';
				break;
		}
		e >>= 2;
	}
	fprintf(stderr, "s%", string);
}
snp_kmer* GetPosFromDictionarySNP(uint8_t isAmbiguos, uint32_t pos, uint8_t altSnpPosInsideKmer, snp_aux_table* snp_aux_tab)
{
	snp_kmer* sn = (snp_kmer*) malloc(sizeof(snp_kmer));

	if (isAmbiguos == 0)
	{
		//fprintf(stderr, "Singolo->NO AUXY TAB.\n");
		if (pos <= 0)
		{
			//fprintf(stderr, "Error: position inside reference genome can't be <= 0");
			exit(0);
		}
		//fprintf(stderr, "Singolo->NO AUXY TAB.\n");
		//fprintf(stderr, "La posizione dentro è: %" PRIu32 "\n", pos);
		sn->positions_to_ref = (uint32_t*) malloc(sizeof(uint32_t));
		sn->posAlternativeSnp = (uint8_t*) malloc(sizeof(uint8_t));

		sn->positions_to_ref[0] = pos;
		sn->posAlternativeSnp[0] = altSnpPosInsideKmer;

		sn->num_occur = 1;
		return sn;
	}
	else if (isAmbiguos == 1)		//REMEMBER: pos_ref_genome is the index of the auxiliary table
	{
		//fprintf(stderr, "Si AUXY TAB\n");
		//fprintf(stderr, "index to acces auxiliary table  %d\n", pos);
		uint8_t i = 0;
		sn->positions_to_ref = (uint32_t*) malloc(10 * sizeof(uint32_t));
		sn->posAlternativeSnp = (uint8_t*) malloc(10 * sizeof(uint8_t));

		//fprintf(stderr,"pos to acces aux table: %"PRIu32 "\n", pos);

		while (true)
		{
			uint32_t currentPos = snp_aux_tab[pos].pos_list[i];
			uint8_t altPosSnpInsideKmer = snp_aux_tab[pos].snp_list[i];
			//i++

			//fprintf(stderr,"current pos read: %"PRIu32 "\n", currentPos);
			//fprintf(stderr,"altPosSnpInsideKmer: %"PRIu8 "\n", altPosSnpInsideKmer);

			if (currentPos == 0)
			{
				//All the positions are been reads
				break;
			}
			else
			{
				sn->positions_to_ref[i] = currentPos;
				sn->posAlternativeSnp[i] = altPosSnpInsideKmer;
				i++;
				if (i == 10)
					break;
			}
		}

		sn->num_occur = i;
		//fprintf(stderr, "Occurence returns %" PRIu8 "\n", i);
		return sn;
	}
	else
	{
		fprintf(stderr, "isAmbiguos assume this impossible value: %" PRIu8 "\n", isAmbiguos);
		exit(0);
	}
}
