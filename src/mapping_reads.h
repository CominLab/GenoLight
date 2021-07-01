#include <cstdint>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include "loadSmartDictionary.h"
#include "utility.h"
#include "fq.hpp"
#include "utility.h"
#include "genolight.h"

#define INDEX_TABLE_SLOT_COUNT  1009
#define INDEX_TABLE_ENTRY_DEPTH  500  /* enough for 5 k-mers */
#define NO_MODIFICATION 10086
#define MAX_HITS 2000

#define ERR_NUM_BASE_MATE_P 25
#define INCREASE_BONUS_MATE_P 1

typedef std::pair<std::size_t, std::size_t> chunck_view;

/*
 * Given a position within the reference genome, I indicate with freq the number
 * of reads that support it (i.e. that the best initial position of that reads is proper index)
 */
typedef struct
{
		uint32_t index;
		uint8_t freq;
} IndexTableEntry;

typedef struct
{
		unsigned short count = 0;
		IndexTableEntry entries[INDEX_TABLE_ENTRY_DEPTH];
} IndexTableSlot;

typedef struct
{
		IndexTableEntry *best;  // highest frequency
		bool ambiguous;         // `best` ambiguous?
		IndexTableSlot table[INDEX_TABLE_SLOT_COUNT];
} IndexTable;

/* convenient way to store k-mer information */
typedef struct
{
		kmer_t kmer;
		uint32_t position;  // 1-based position of read based on kmer hit
		uint32_t kmer_pos;  // 1-based position of k-mer
		uint32_t modified_pos; // 0-based position in range of [0.32)
		uint8_t altSnpPosInsideKmer;
		bool is_neighbor;

} kmer_context;

typedef struct
{
		//Position inside pileup table previously updated
		uint32_t posPileupTable;
		/*
		 * This variable is true if the updated base is the reference one,
		 * false if it is the alternative one
		 */
		bool isRef;

} infoUpdatePileupTable;


void updatePilleupSingleEnd(param_creation_geno* geno_param, def_param* dp);
void updatePilleupMatePairs(param_creation_geno* geno_param, def_param* dp);
