#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include "pileup.h"
#include <inttypes.h>
#include <dictgen.h>
#include <math.h>
#define LOAD_FACTOR 0.4

double near_power_of_two(const uint32_t number);

TotaleKmerSnp* total_kmer_snp_struct_init(const size_t numerKmer)
{
	TotaleKmerSnp* p = (TotaleKmerSnp*) malloc(sizeof(TotaleKmerSnp));
	const size_t size = pow(2, near_power_of_two(numerKmer));
	assert((size & (size - 1)) == 0);  // `size` must be a power of 2
	fprintf(stderr, "Dimensione total kmer snp is: %zu \n", size);
	p->count = 0;
	p->size = size;
	p->threshold = (size_t)(size * LOAD_FACTOR);

	p->table = (struct total_kmer_snp_entry**) calloc(size, sizeof(*(p->table)));
	assert(p->table);
	return p;
}

void total_kmer_snp_struct_dealloc(TotaleKmerSnp *p)
{
	struct total_kmer_snp_entry **table = p->table;
	const size_t size = p->size;
	for (size_t i = 0; i < size; i++)
	{
		struct total_kmer_snp_entry *e = table[i];
		while (e != NULL)
		{
			struct total_kmer_snp_entry *temp = e;
			e = e->next;
			free(temp->val);
			free(temp);
		}
	}
	free(table);
}

void total_kmer_snp_struct_add(TotaleKmerSnp *p, const uint32_t key, const uint32_t pos, const uint8_t snp)
{
	value* v = total_kmer_snp_struct_get(p, key);

	if (v != NULL)
	{
		v->pos = pos;
		v->snpInfo = snp;
		return;
	}
	const size_t size = p->size;
	struct total_kmer_snp_entry **table = p->table;
	const uint32_t n = hashFunction(key) & (size - 1);

	struct total_kmer_snp_entry* e = (struct total_kmer_snp_entry *) malloc(sizeof(*e));
	//e->pos = pos;
	//e->key = key;

	v = (value*) malloc(sizeof(value));
	v->pos = pos;
	v->snpInfo = snp;
	e->val = v;
	e->key = key;
	e->next = table[n];
	table[n] = e;
	/*
	 if (p->count++ > p->threshold)
	 {
	 grow(p);
	 }
	 */
}

PrincipalDictionaty* principalDictionary_init(const uint32_t sizeDictionaty)
{
	PrincipalDictionaty* p = (PrincipalDictionaty*) malloc(sizeof(PrincipalDictionaty));
	uint32_t size = pow(2, near_power_of_two(sizeDictionaty));

	fprintf(stderr, "Dimensione dizionario principale è: %zu \n", size);

	assert((size & (size - 1)) == 0);  // `size` must be a power of 2

	p->count = 0;
	p->size = size;
	//p->threshold = (size_t)(size * LOAD_FACTOR);

	p->table = (struct principal_entry**) calloc(size, sizeof(*(p->table)));
	assert(p->table);
	return p;
}

void PrincipalDictionaradd(PrincipalDictionaty *p, const uint32_t key, const value_snp* v)
{
	const uint32_t n = hashFunction(key) & (p->size - 1);

	struct principal_entry *e = p->table[n];

	//Table don't have an element inside position n, so i create the first one
	if (e == NULL)
	{
		e = (principal_entry*) malloc(sizeof(principal_entry));
		e->value = v;
		e->key = key;
		e->next = NULL;
		p->table[n] = e;
	}
	else
	{
		//Table have an element inside position n,so i go to the end of bucket and
		//store
		while (true)
		{
			if (e->next != NULL)
				e = e->next;
			else
			{
				struct principal_entry *m = (principal_entry*) malloc(sizeof(principal_entry));
				m->value = v;
				m->key = key;
				m->next = NULL;

				e->next = m;
				break;
			}
		}
	}
}

void principalDictionaty_dealloc(PrincipalDictionaty *p)
{
	struct principal_entry **table = p->table;

	const size_t size = p->size;

	for (size_t i = 0; i < size; i++)
	{
		struct principal_entry *e = table[i];

		while (e != NULL)
		{
			free(e->value);
			struct principal_entry *temp = e;
			e = e->next;
			free(temp);
		}
	}
	free(table);
}

/*
 * Method that ritorna un numero inferiore alla potenza di 2 più piccolo
 */
double near_power_of_two(const uint32_t number)
{
	double i = 0;

	//Only for consistency..
	if (number <= 2)
		return 2;

	while (true)
	{
		i = i + 1;
		double p = pow(2, i);

		if ((number - p) < 0)
			return i - 2;
	}
}

void ptable_init(PileupTable *p, const size_t size)
{
	assert((size & (size - 1)) == 0);  // `size` must be a power of 2
	p->table = (struct pileup_entry**) calloc(size, sizeof(*(p->table)));
	assert(p->table);
	p->count = 0;
	p->size = size;
	p->threshold = (size_t)(size * LOAD_FACTOR);
}

void ptable_dealloc(PileupTable *p)
{
	struct pileup_entry **table = p->table;
	const size_t size = p->size;
	for (size_t i = 0; i < size; i++)
	{
		struct pileup_entry *e = table[i];
		while (e != NULL)
		{
			struct pileup_entry *temp = e;
			e = e->next;
			free(temp);
		}
	}
	free(table);
}

static void grow(PileupTable *p)
{
	const size_t size = p->size;
	struct pileup_entry **table = p->table;

	const size_t new_size = 2 * size;
	struct pileup_entry **new_table = (struct pileup_entry**) calloc(new_size, sizeof(*new_table));
	assert(new_table);

	for (size_t i = 0; i < size; i++)
	{
		struct pileup_entry *e = table[i];
		while (e != NULL)
		{
			struct pileup_entry *next = e->next;
			const uint32_t n = hashFunction(e->key) & (new_size - 1);
			e->next = new_table[n];
			new_table[n] = e;
			e = next;
		}
	}

	p->table = new_table;
	p->size = new_size;
	p->threshold = (size_t)(new_size * LOAD_FACTOR);
	free(table);
}

void ptable_add(PileupTable *p, const uint32_t key, unsigned ref, unsigned alt, uint8_t ref_freq, uint8_t alt_freq)
{
	if (ptable_get(p, key) != NULL)
		return;

	const size_t size = p->size;
	struct pileup_entry **table = p->table;

	const uint32_t n = hashFunction(key) & (size - 1);

	struct pileup_entry *e = (struct pileup_entry *) malloc(sizeof(*e));
	e->ref = ref;
	e->alt = alt;
	e->ref_cnt = 0;
	e->alt_cnt = 0;
	e->ref_freq = ref_freq;
	e->alt_freq = alt_freq;
	e->key = key;
	e->next = table[n];
	table[n] = e;

	if (p->count++ > p->threshold)
	{
		grow(p);
	}
}

