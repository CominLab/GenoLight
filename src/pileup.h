#ifndef PILEUP_H
#define PILEUP_H

#include <stdlib.h>
#include <stdint.h>

struct snp_aux_table
{
		uint32_t* pos_list;
		uint8_t* snp_list;
}__attribute__((packed));

struct Entry_dictionary_snp
{
		uint32_t pos_snp_genome = 0;
		uint8_t is_reverse :1;
		uint8_t is_compleate :1;
		uint8_t isDirectAmb :1;
		uint8_t isRevAmb :1;
		uint8_t altSnpPosInsideKmer;
		uint8_t alt_base :3;
}__attribute__((packed));

struct value_snp
{
		uint32_t pos_snp_genome = 0;
		uint8_t is_reverse :1;
		uint8_t is_compleate :1;
		uint8_t isDirectAmb :1;
		uint8_t isRevAmb :1;
		uint8_t altSnpPosInsideKmer;
		uint8_t alt_base :3;
}__attribute__((packed));

//ENTRY UNORDERED MAP
struct principal_entry
{
		//Value
		value_snp* value;
		//Key
		uint32_t key;
		//Pointer to the next reference
		struct principal_entry *next;
};

//UNORDERED MAP DIZIONARIO PRINCIPALE
typedef struct
{
		struct principal_entry **table;
		size_t count;
		size_t size;
		size_t threshold;
} PrincipalDictionaty;

//------------------------------------------------------------------

//Struct that contains information store in unordered map
struct value
{
		//Posizione all'intero del genoma modificato contenente	la k-mer r.c. di una k-mer completa
		uint32_t pos;
		//Contiene la posizione dello snp all'interno di tale k-mer
		uint8_t snpInfo;

}__attribute__((packed));

struct total_kmer_snp_entry
{
		struct value* val;
		uint32_t key;
		struct total_kmer_snp_entry *next;
};

typedef struct
{
		struct total_kmer_snp_entry **table;
		size_t count;
		size_t size;
		size_t threshold;
} TotaleKmerSnp;

struct pileup_entry
{
		unsigned ref :2;
		unsigned alt :2;
		unsigned ref_cnt :6;
		unsigned alt_cnt :6;
		uint8_t ref_freq;
		uint8_t alt_freq;
		uint32_t key;
		struct pileup_entry *next;
};

typedef struct
{
		struct pileup_entry **table;
		size_t count;
		size_t size;
		size_t threshold;
} PileupTable;

/*
 * Analogamente anche il dizionario degli SNP è composto da 3 componenti che sono:
 */
typedef struct
{
		//Principale
		//PrincipalDictionaty* snp_dict;
		Entry_dictionary_snp* snp_dict;
		uint32_t sizeDictionary;

		//Contiene le info dei duplicati
		struct snp_aux_table* snp_auxiliary_table;
		uint32_t sizeAuxiliaryTable;

		TotaleKmerSnp* reverse_total_kmer_snp;

} SNPSmartDictionary;

TotaleKmerSnp* total_kmer_snp_struct_init(const size_t size);
PrincipalDictionaty* principalDictionary_init(const uint32_t sizeDictionary);

void total_kmer_snp_struct_dealloc(TotaleKmerSnp *p);

void total_kmer_snp_struct_add(TotaleKmerSnp *p, const uint32_t key, const uint32_t pos, const uint8_t snpInfo);
void PrincipalDictionaradd(PrincipalDictionaty *p, const uint32_t key, const value_snp* v);

uint32_t* GetPosFromDictionarySNP(uint8_t isAmbiguose, uint32_t pos, snp_aux_table* snp_aux_tab, uint8_t* real_num_occ);
void principalDictionaty_dealloc(PrincipalDictionaty *p);

void ptable_init(PileupTable *p, const size_t size);
void ptable_dealloc(PileupTable *p);
void ptable_add(PileupTable *p, const uint32_t key, unsigned ref, unsigned alt, uint8_t ref_freq, uint8_t alt_freq);

/*
 * Adapted from java.util.HashMap
 */
static inline uint32_t hashFunction(uint32_t h)
{
	h ^= (h >> 20) ^ (h >> 12);
	return h ^ (h >> 7) ^ (h >> 4);
}

static inline value* total_kmer_snp_struct_get(TotaleKmerSnp *p, const uint32_t key)
{
	for (struct total_kmer_snp_entry *e = p->table[hashFunction(key) & (p->size - 1)]; e != NULL; e = e->next)
	{
		if (e->key == key)
		{
			return e->val;
		}
	}
	return NULL;
}

static inline struct value_snp* principalDictionary_get(PrincipalDictionaty *p, const uint32_t key)
{
	for (struct principal_entry *e = p->table[hashFunction(key) & (p->size - 1)]; e != NULL; e = e->next)
	{
		if (e->key == key)
		{
			return e->value;
		}
	}
	return NULL;
}

static inline struct pileup_entry *ptable_get(PileupTable *p, const uint32_t key)
{
	for (struct pileup_entry *e = p->table[hashFunction(key) & (p->size - 1)]; e != NULL; e = e->next)
	{
		if (e->key == key)
		{
			return e;
		}
	}
	return NULL;
}

#endif /* PILEUP_H */

