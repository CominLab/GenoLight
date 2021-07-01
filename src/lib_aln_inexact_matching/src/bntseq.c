/* The MIT License

 Copyright (c) 2008 Genome Research Ltd (GRL).

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
#include <zlib.h>
#include <unistd.h>
#include <errno.h>
#include <stdbool.h>
#include "bntseq.h"
#include "utils.h"

#include "kseq.h"
KSEQ_DECLARE(gzFile)

#include "khash.h"
KHASH_MAP_INIT_STR(str, int)

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

unsigned char nst_nt4_table[256] =
{4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5 /*'-'*/, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4};


static uint8_t *add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q);
static void bns_dump(const bntseq_t *bns, const char *prefix);
static void bns_destroy(bntseq_t *bns);

/*
 * Derived version of bwa_sa2pos present in bwase.c
 */
uint64_t bwa_sa2pos(const bntseq_t *bns, const bwt_t *bwt, uint64_t sapos, int ref_len, bool *strand)
{
	uint64_t pos_f;
	bool is_rev;
	*strand = 0; // initialise strand to 0 otherwise we could return without setting it
	pos_f = bwt_sa(bwt, sapos); // position on the forward-reverse coordinate
	if (pos_f < bns->l_pac && bns->l_pac < pos_f + ref_len)
		return (uint64_t) -1;
	pos_f = bns_depos(bns, pos_f, &is_rev); // position on the forward strand; this may be the first base or the last base
	*strand = !is_rev;
	if (is_rev)
		pos_f = pos_f + 1 < ref_len ? 0 : pos_f - ref_len + 1; // position of the first base
	return pos_f + 1; // FIXME: it is possible that pos_f < bns->anns[ref_id].offset
}

uint8_t *bns_get_seq(int64_t l_pac, const uint8_t *pac, int64_t beg, int64_t end, int64_t *len)
{
	beg--;//base 0
	end--;
	uint8_t *seq = 0;
	if (end < beg)
		end ^= beg, beg ^= end, end ^= beg; // if end is smaller, swap
	if (end > l_pac << 1)
		end = l_pac << 1;
	if (beg < 0)
		beg = 0;

	if (beg >= l_pac || end <= l_pac)
	{
		int64_t k, l = 0;
		*len = end - beg;
		seq = malloc(end - beg);
		if (beg >= l_pac)
		{ // reverse strand
			int64_t beg_f = (l_pac << 1) - 1 - end;
			int64_t end_f = (l_pac << 1) - 1 - beg;
			for (k = end_f; k > beg_f; --k)
				seq[l++] = 3 - _get_pac(pac, k);
		}
		else
		{ // forward strand
			for (k = beg; k < end; ++k)
				seq[l++] = _get_pac(pac, k);
		}
	}
	else
		*len = 0; // if bridging the forward-reverse boundary, return nothing
	return seq;
}

int64_t bns_fasta2bntseq(gzFile fp_fa, const char *prefix, int for_only)
{
	extern void seq_reverse(int len, ubyte_t *seq, int is_comp); // in bwaseqio.c
	kseq_t *seq;
	char name[1024];
	bntseq_t *bns;
	uint8_t *pac = 0;
	int32_t m_seqs, m_holes;
	int64_t ret = -1, m_pac, l;
	bntamb1_t *q;
	FILE *fp;

	// initialization
	seq = kseq_init(fp_fa);
	bns = (bntseq_t*) calloc(1, sizeof(bntseq_t));
	bns->seed = 11; // fixed seed for random generator
	srand48(bns->seed);
	m_seqs = m_holes = 8;
	m_pac = 0x10000;
	bns->anns = (bntann1_t*) calloc(m_seqs, sizeof(bntann1_t));
	bns->ambs = (bntamb1_t*) calloc(m_holes, sizeof(bntamb1_t));
	pac = calloc(m_pac / 4, 1);
	q = bns->ambs;
	strcpy(name, prefix);
	strcat(name, ".pac");
	fp = xopen(name, "wb");
	// read sequences
	while (kseq_read(seq) >= 0)
		pac = add1(seq, bns, pac, &m_pac, &m_seqs, &m_holes, &q);
	if (!for_only)
	{ // add the reverse complemented sequence
		int64_t ll_pac = (bns->l_pac * 2 + 3) / 4 * 4;
		if (ll_pac > m_pac)
			pac = realloc(pac, ll_pac / 4);
		memset(pac + (bns->l_pac + 3) / 4, 0, (ll_pac - (bns->l_pac + 3) / 4 * 4) / 4);
		for (l = bns->l_pac - 1; l >= 0; --l, ++bns->l_pac)
			_set_pac(pac, bns->l_pac, 3-_get_pac(pac, l));
	}
	ret = bns->l_pac;
	{ // finalize .pac file
		ubyte_t ct;
		err_fwrite(pac, 1, (bns->l_pac >> 2) + ((bns->l_pac & 3) == 0 ? 0 : 1), fp);
		// the following codes make the pac file size always (l_pac/4+1+1)
		if (bns->l_pac % 4 == 0)
		{
			ct = 0;
			err_fwrite(&ct, 1, 1, fp);
		}
		ct = bns->l_pac % 4;
		err_fwrite(&ct, 1, 1, fp);
		// close .pac file
		err_fflush(fp);
		err_fclose(fp);
	}
	bns_dump(bns, prefix);
	bns_destroy(bns);
	kseq_destroy(seq);
	free(pac);
	return ret;
}

static uint8_t *add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q)
{
	bntann1_t *p;
	int i, lasts;
	if (bns->n_seqs == *m_seqs)
	{
		*m_seqs <<= 1;
		bns->anns = (bntann1_t*) realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
	}
	p = bns->anns + bns->n_seqs;
	p->name = strdup((char*) seq->name.s);
	p->anno = seq->comment.l > 0 ? strdup((char*) seq->comment.s) : strdup("(null)");
	p->gi = 0;
	p->len = seq->seq.l;
	p->offset = (bns->n_seqs == 0) ? 0 : (p - 1)->offset + (p - 1)->len;
	p->n_ambs = 0;
	for (i = lasts = 0; i < seq->seq.l; ++i)
	{
		int c = nst_nt4_table[(int) seq->seq.s[i]];
		if (c >= 4)
		{ // N
			if (lasts == seq->seq.s[i])
			{ // contiguous N
				++(*q)->len;
			}
			else
			{
				if (bns->n_holes == *m_holes)
				{
					(*m_holes) <<= 1;
					bns->ambs = (bntamb1_t*) realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
				}
				*q = bns->ambs + bns->n_holes;
				(*q)->len = 1;
				(*q)->offset = p->offset + i;
				(*q)->amb = seq->seq.s[i];
				++p->n_ambs;
				++bns->n_holes;
			}
		}
		lasts = seq->seq.s[i];
		{ // fill buffer
			if (c >= 4)
				c = lrand48() & 3;
			if (bns->l_pac == *m_pac)
			{ // double the pac size
				*m_pac <<= 1;
				pac = realloc(pac, *m_pac / 4);
				memset(pac + bns->l_pac / 4, 0, (*m_pac - bns->l_pac) / 4);
			}
			_set_pac(pac, bns->l_pac, c);
			++bns->l_pac;
		}
	}
	++bns->n_seqs;
	return pac;
}

static void bns_dump(const bntseq_t *bns, const char *prefix)
{
	char str[1024];
	FILE *fp;
	int i;
	{ // dump .ann
		strcpy(str, prefix);
		strcat(str, ".ann");
		fp = xopen(str, "w");
		err_fprintf(fp, "%lld %d %u\n", (long long) bns->l_pac, bns->n_seqs, bns->seed);
		for (i = 0; i != bns->n_seqs; ++i)
		{
			bntann1_t *p = bns->anns + i;
			err_fprintf(fp, "%d %s", p->gi, p->name);
			if (p->anno[0])
				err_fprintf(fp, " %s\n", p->anno);
			else
				err_fprintf(fp, "\n");
			err_fprintf(fp, "%lld %d %d\n", (long long) p->offset, p->len, p->n_ambs);
		}
		err_fflush(fp);
		err_fclose(fp);
	}
	{ // dump .amb
		strcpy(str, prefix);
		strcat(str, ".amb");
		fp = xopen(str, "w");
		err_fprintf(fp, "%lld %d %u\n", (long long) bns->l_pac, bns->n_seqs, bns->n_holes);
		for (i = 0; i != bns->n_holes; ++i)
		{
			bntamb1_t *p = bns->ambs + i;
			err_fprintf(fp, "%lld %d %c\n", (long long) p->offset, p->len, p->amb);
		}
		err_fflush(fp);
		err_fclose(fp);
	}
}

static void bns_destroy(bntseq_t *bns)
{
	if (bns == 0)
		return;
	else
	{
		int i;
		if (bns->fp_pac)
			err_fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i)
		{
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}
