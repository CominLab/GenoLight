#include <string>
#include "lib_aln_inexact_matching.h"
#include "mapping_reads.h"
#include "fq.hpp"
#include "geno.h"
#include "utility.h"
#include "CTPL/ctpl_stl.h"
#include "semaphore.h"
#include "inttypes.h"
#include <unistd.h>
#include <chrono>
#include<cstdlib>

typedef std::pair<std::size_t, std::size_t> chunck_view;

const int verbose_map = 0;

double startLib, endLib, timeMedio = 0;
int readN = 0;

static void UpdateIndexTable(IndexTable* index_table, search_result* kmer_neighborhood, kmer_context* hit_contexts,
								unordered_map<uint32_t, unordered_set<uint32_t>> & index_2_kmer_pos_set, size_t* n_hi,
								param_creation_geno* geno_param);

static void improved_index_table_add(IndexTable *index_table, uint32_t index, uint32_t kmer_pos,
										unordered_map<uint32_t, unordered_set<uint32_t>> & index_2_kmer_pos_set, bool is_neighbor);

static void index_table_clear_index(IndexTable *index_table, uint32_t index);

static void index_table_clear(IndexTable *index_table);

static void UpdatePileupTable(kmer_context* hit_contexts, size_t n_hits, uint32_t target_index, param_creation_geno* geno_param,
								std::vector<infoUpdatePileupTable>& infoVect, IndexTable* index_table);

void partialUpdate(std::vector<infoUpdatePileupTable> snpToUpdate, param_creation_geno* geno_param);

void UpdateSingleEndReads(param_creation_geno* geno_param, std::vector<FastqRecord> &chunk, def_param* dp);

/**
 * Core method in which attempts to map the read in the reference genome and, if it fails, tries to align the r.c of
 * the reads. In reality a real alignment is not carried out, but the position within the reference genome supported by
 * several k-mer extracted from the reads is searched. This is done using an IndexTable as the data structure.
 */
void UpdateSingleEndReads(param_creation_geno* geno_param, std::vector<FastqRecord> &chunk, def_param* dp)
{
	std::vector<infoUpdatePileupTable> updatePileupTable;
	output* search_output = nullptr;

	/**
	 * Data structure used to identify the position supported by current read
	 * (that is the most probable position in which the read aligns whit the
	 * reference genome)
	 */
	//std::cerr << "Index table ..." << std::endl;
	IndexTable index_table;
	index_table_clear(&index_table);

	/**
	 * It tracks the position of each k-mer within the reference genome(obviously it is reset at each new reads).
	 */
	kmer_context hit_contexts[MAX_HITS * 2];
	size_t n_hits = 0;
	//std::cerr << "for ..." << std::endl;
	//For each record(a record contains sequence of bases and quality score of a read) contained in the single chunk
	for (int e = 0; e < chunk.size(); e++)
	{
		//I recover a record of thread's chunk
		//std::cerr << "1..." << std::endl;
		FastqRecord record = chunk[e];

		n_hits = 0;
		//std::cerr << "2..." << std::endl;
		unordered_map<uint32_t, unordered_set<uint32_t>> index_2_kmer_pos_set;
		readN++;
		//std::cerr << "3..." << std::endl;
		std::string rec = record.seq;
		//std::cerr << "4 ..." << std::endl;
		//Number of hits found by lib_aln_bound_backtracking
		search_output = lib_aln_bound_backtracking(dp, record.seq.c_str(), record.qua.c_str());
		//std::cerr << "if ..." << std::endl;
		//No hit found
		if (search_output->n_for + search_output->n_rev != 0)
		{

			//I'm sure at least one hit must be there
			if (verbose_map > 5)
				std::cerr << "\nThe number of hit found is: " << search_output->n_for + search_output->n_rev << std::endl;

			//Iteratively consider all the k-mer present in result_ref
			for (uint32_t l = 0; l < search_output->n_for; l++)
			{
				if (search_output->res_for[l]->num_occur > 10)
					continue;

				if (verbose_map > 100)
				{
					std::cerr << "The current k-mer DIRECT-REFERENCE (return by library) analyzed is: "
					<< decode_kmer_vargeno(search_output->res_for[l]->hit, 32) << std::endl;

					std::cerr << "---->num_occ: " << (unsigned) search_output->res_for[l]->num_occur << std::endl;

					for (int i = 0; i < search_output->res_for[l]->num_occur; i++)
					{
						std::cerr << "Position_" << i << "- " << search_output->res_for[l]->positions_to_ref[i] << std::endl;
					}
					std::cerr << "offset " << (unsigned) search_output->res_for[l]->offset << std::endl;

					std::cerr << "n_mismatch: " << unsigned(search_output->res_for[l]->n_mismatches) << std::endl;

					if (search_output->res_for[l]->n_mismatches != 0)
						std::cerr << "different_positions: " << unsigned(search_output->res_for[l]->different_positions[0]) << std::endl;
				}

				UpdateIndexTable(&index_table, search_output->res_for[l], hit_contexts, index_2_kmer_pos_set, &n_hits, geno_param);
			}

			for (uint32_t l = 0; l < search_output->n_rev; l++)
			{
				if (search_output->res_rev[l]->num_occur > 10)
					continue;

				if (verbose_map > 100)
				{
					std::cerr << "The current k-mer REVERSE (return by library) analyzed is: "
					<< decode_kmer_vargeno(search_output->res_rev[l]->hit, 32) << std::endl;

					std::cerr << "---->num_occ: " << (unsigned) search_output->res_rev[l]->num_occur << std::endl;

					for (int i = 0; i < search_output->res_rev[l]->num_occur; i++)
					{
						std::cerr << "Position_" << i << "- " << search_output->res_rev[l]->positions_to_ref[i] << std::endl;
					}

					std::cerr << "offset " << (unsigned) search_output->res_rev[l]->offset << std::endl;

					std::cerr << "n_mismatch: " << unsigned(search_output->res_rev[l]->n_mismatches) << std::endl;

					if (search_output->res_rev[l]->n_mismatches != 0)
						std::cerr << "different_positions: " << unsigned(search_output->res_rev[l]->different_positions[0]) << std::endl;
				}
			}
		}

		bool process_read = index_table.best && (index_table.best->freq > 1) && !index_table.ambiguous;

		if (process_read == false)
		{
			//In this code block 'process_read' can become true if we map the r.c. of the read
			index_table.best = NULL;
			index_table.ambiguous = false;

			n_hits = 0;

			for (int j = search_output->n_rev - 1; j >= 0; j--)
			{
				if (search_output->res_rev[j]->num_occur > 10)
					continue;

				//Scandisco ogni hit trovata e aggiorno l'index table
				UpdateIndexTable(&index_table, search_output->res_rev[j], hit_contexts, index_2_kmer_pos_set, &n_hits, geno_param);
			}

			process_read = (index_table.best && (index_table.best->freq > 1) && !index_table.ambiguous);
		}

		if (process_read == true)
		{
			//std::cerr << "****(308)La read " << record->seq << " è presente (diretta o r.c)****" << std::endl;

			//If there is a best, then I take the best index (now the value of freq field I do not care anymore)
			const uint32_t target_index = index_table.best ? index_table.best->index : 0;

			if (verbose_map > 2)
				std::cerr << "\n\n***UpdatePileupTable***" << std::endl;

			if (target_index != 0)
			{
				if (verbose_map > 2)
				{
					std::cerr << "n_hits: " << n_hits << std::endl;
					std::cerr << "\ntarget_index: " << target_index << std::endl;
				}
				UpdatePileupTable(hit_contexts, n_hits, target_index, geno_param, updatePileupTable, &index_table);
			}
		}
		else
		{
			if (verbose_map > 2)
				std::cerr << "****(323)La read " << record.seq << " NON è presente (diretta o r.c)****" << std::endl;
		}

		lib_aln_sr_destroy(search_output);
		search_output = nullptr;

		//In this code block 'process_read' can become true if we map the r.c. of the read
		index_table.best = NULL;
		index_table.ambiguous = false;

		//Pulisco l'index table
		for (size_t i = 0; i < n_hits; i++)
		{
			const uint32_t index = hit_contexts[i].position;
			index_table_clear_index(&index_table, index);
		}

		//std::cerr << "Terminato mappare una reads" << std::endl;

	}	//All the reads have been used
}


/**
 *@geno_param: reference to a struct that contains all information concerning the genotyping phase
 */
void updatePilleupSingleEnd(param_creation_geno* geno_param, def_param* dp)
{
	//Get chunck size specified by user. Default is 1000
	int chunk_size = geno_param->chunk_size;
	std::cerr << "chunk_size: " << chunk_size << std::endl;

	//Counter
	int numerReadsAnalyze = 0;

	//Counter of number of reads analyzed by program
	uint64_t totalReadAnalyzed = 0;

	//Time information
	clock_t partial_begin, partial_end;

	//Open a flow to file contains reads by creating an object of type fq named read_fastq
	fq* read_fastq = new fq(geno_param->path_fastq);

	//Circular vector of size two contains two vector. This two vector stored one the current
	//chunck analyzed by the program and the second one the chunck that will be process
	//next interation
	std::vector <FastqRecord> chunks(chunk_size);

	//Set which chunck is elaborated. Can be see as the pointer to the circular vector.
	uint8_t chunk_counter = 0;

	const clock_t begin = clock();

	partial_begin = begin;

	std::cerr << "Start time" << debugProgress() << std::endl;

	do
	{
		//std::cerr << "Fill chunk start!" << std::endl;

		//Vector containing geno_param-> chunk size reads
		read_fastq->fill_chunk(chunks);
		//std::cerr << "Fill chunk end!" << std::endl;

		if (chunks.size() != 0)
		{
			//std::cerr << "Update Single End Reads" << std::endl;
			UpdateSingleEndReads(geno_param, std::ref(chunks), dp);
			totalReadAnalyzed += chunks.size();
		}
		else
			std::cerr << "chunk vuoto" << std::endl;

	}
	while (read_fastq->AllReadExtracted() == false);

	delete read_fastq;

	clock_t const end = clock();
	double const time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

	std::cout << "\nNumero di reads analizzate: " << totalReadAnalyzed << std::endl;

	printf("(446)-Tempo per effettuare il mapping: %f sec\n", time_spent);
	std::cout << "End mapping reads!\n" << std::endl;
}


/**
 *@index_table
 *@kmer_neighborhood
 *@hit_contexts
 *@index_2_kmer_pos_set
 *@pileup_table
 *@pileup_table
 *@n_hi
 *@geno_param: reference to a struct that contains all information concerning the genotyping phase
 */
static void UpdateIndexTable(IndexTable* index_table, search_result* kmer_neighborhood, kmer_context* hit_contexts,
								unordered_map<uint32_t, unordered_set<uint32_t>> & index_2_kmer_pos_set, size_t *n_hi,
								param_creation_geno* geno_param)
{
	//std::cerr << "\n*****(475)-UpdateIndexTable" << std::endl;

	int n_hits = *(n_hi);

	//std::vector < uint32_t > position = kmer_neighborhood->pos_genome;
	if ((kmer_neighborhood->num_occur == 1) && (kmer_neighborhood->positions_to_ref[0] == 0))
	{
		return;
	}
	/*
	 * Iteratively consider all the positions in which 'kmer_neighborhood.kmer' appears within the reference genome.
	 * In most cases, k-mer is present only once within the reference genome and therefore the vector
	 * pos_kmer_genome has only one element.
	 */
	for (std::vector<unsigned int>::size_type j = 0; j < kmer_neighborhood->num_occur; j++)
	{
		//I deduce the position of the corresponding reads to which the k-mer belongs
		const uint32_t read_pos = kmer_neighborhood->positions_to_ref[j] - kmer_neighborhood->offset;

		if (kmer_neighborhood->n_mismatches == 0)
		{
			//Current k-mer is equal to kmer_to_search. For this reason NO_MODIFICATION
			kmer_context k_temp
			{kmer_neighborhood->hit, read_pos, kmer_neighborhood->positions_to_ref[j], NO_MODIFICATION};

			hit_contexts[n_hits++] = k_temp;

			if (verbose_map > 1)
				std::cerr << "This k-mer is present without mismatch inside reference genome" << std::endl;
			improved_index_table_add(index_table, read_pos, kmer_neighborhood->positions_to_ref[j], index_2_kmer_pos_set, false);
		}
		else
		{
			//The current k-mer has a hamming distance of 1 compared to 'kmer_to_search'..

			//kmer_neighborhood.pos_mismatch is 1-base but in this case i need 0-base
			const uint32_t diff_base_pos = kmer_neighborhood->different_positions[0] - 1;

			if ((diff_base_pos < 0) || (diff_base_pos >= 32))
			{
				std::cerr << "Diff pos must be between 0-31!" << std::endl;
				std::cerr << "Diff-base-pos: " << (unsigned) diff_base_pos << std::endl;
				exit(0);
			}
			//If diff_base_pos is 0-base, then ref_hit_diff_loc is the position inside reference genome
			const size_t ref_hit_diff_loc = kmer_neighborhood->positions_to_ref[j] + diff_base_pos;

			//Se in tale posizione non ci sono SNP
			/*
			 if ((kmer_neighborhood->altSnpPosInsideKmer == NULL)
			 && ((geno_param->pileup_table[ref_hit_diff_loc].ref != 0) || (geno_param->pileup_table[ref_hit_diff_loc].alt != 0)))
			 continue;
			 */

			if ((kmer_neighborhood->altSnpPosInsideKmer == NULL) && (ptable_get(&geno_param->ptable, ref_hit_diff_loc) == NULL))
				continue;

			//COMENTO
			if ((kmer_neighborhood->altSnpPosInsideKmer != NULL) && (kmer_neighborhood->altSnpPosInsideKmer[j] == diff_base_pos))
				continue;

			kmer_context k_temp
			{kmer_neighborhood->hit, read_pos, kmer_neighborhood->positions_to_ref[j], diff_base_pos};

			if (verbose_map > 1)
			{
				std::cerr << "\nRef_hit_diff_loc: " << unsigned(ref_hit_diff_loc) << std::endl;
				std::cerr << "Diff_base_pos: " << unsigned(diff_base_pos) << std::endl;
			}

			hit_contexts[n_hits++] = k_temp;

			improved_index_table_add(index_table, read_pos, kmer_neighborhood->positions_to_ref[j], index_2_kmer_pos_set, true);

		}
	}
	*n_hi = n_hits;
}

/**
 * Method in which the pileup table is updated.
 *
 * @hit_contexts: Pointer to an array containing the information regarding the k-mer that belong to the current read.
 * 				  As mentioned many times, such k-mer can really be the k-mer (not overlapping) extracted from the reads
 * 				  or belong to N (1) of a non overlapping k-mer extracted from the reads.
 *
 * @n_hits: Number of k-mer inside vector
 *
 * @target_index: Position within the genome in which the reads is positioned
 *
 * @pileup_table: Reference to pileup table
 *
 * @pileup_size: Size of pileup table
 *
 * @k: Size of k-mer
 */
static void UpdatePileupTable(kmer_context* hit_contexts, size_t n_hits, uint32_t target_index, param_creation_geno* geno_param,
								std::vector<infoUpdatePileupTable>& infoVect, IndexTable* index_table)
{
	//std::cerr << "**(586)-UPDATE PILEUP TABLE " << std::endl;
	int numSnpFound = 0;

	//For each k-mer inside vector has pointer hit_contexts
	for (size_t i = 0; i < n_hits; i++)
	{
		//Check if ref_hit_contexts[i] is a k-mer that supports target_index,otherwise I do nothing
		const uint32_t index = hit_contexts[i].position;

		//Pulisco l'index table
		index_table_clear_index(index_table, index);

		if (index == target_index)
		{
			//Start position of k-mer inside the reference genome
			const uint32_t kmer_pos = hit_contexts[i].kmer_pos;

			//Coded sequence of k-mer
			const kmer_t kmer = hit_contexts[i].kmer;

			//std::cerr << "\n--->Current k-mer to update pileup table: " << decode_kmer_vargeno(kmer, geno_param->k) << std::endl;

			//Position in which this k-mer differs with respect to the "original" k-mer present in the reads
			const uint32_t modified_pos = hit_contexts[i].modified_pos;

			/**
			 * Iteratively scan all the bases of current k-mer analyzed.
			 */
			for (unsigned j = 0; j < geno_param->k; j++)
			{
				if ((geno_param->k - j - 1) == modified_pos)
				{
					//std::cerr << "EXIT MODIFIED POS" << (unsigned) (geno_param->k - j - 1) << std::endl;
					continue;
				}
				/**
				 * Example: If k-mer is ACG, at the first iteration method kmer_get_base(ACG,0)
				 * return G so the position within the pileup table that I must consider is
				 * kmer_pos - j + k - 1
				 */
				const unsigned base = kmer_get_base(kmer, j);

				//struct packed_pileup_entry *p = &geno_param->pileup_table[kmer_pos - j + geno_param->k - 1];
				struct pileup_entry* p = ptable_get(&geno_param->ptable, kmer_pos - j + geno_param->k - 1);

				/*
				 * Inside pileup_table some index doesn't match any SNP
				 * (see how pileup table is create in method load_smartSNP_Dict inside
				 * loadSmartDictionary.cc)
				 */
				//std::cerr << "pileup table in position: " << kmer_pos - j + geno_param->k - 1 << std::endl;
				if (p != NULL)
				{
					//This means that a SNP is present in kmer_pos + i
					if (base == p->ref)
					{
						//MAX_COV((1 << 6)-1)
						if (p->ref_cnt != MAX_COV)
							++p->ref_cnt;
					}
					else if (base == p->alt)
					{
						if (p->alt_cnt != MAX_COV)
							++p->alt_cnt;
					}
				}
			}

		}
	}
}

/*
 * Given kmer(parameter), consider all k-mer that belong to N(1). Each one is searched in
 * reassembly(forward and r.c).
 *
 * For each hit, create a struct of type pos_neighborhoodin which store the
 * Attention: the basis corresponding to the positions NAME_DA_DEFINIRE_ will not be modified.
 *
 * NOTE: it's possible that a k-mer isn't found in the REASSEMBLY. It simply means that in the
 * 		 REFERENCE GENOMA there is neither such k-mer and its r.c.
 *
 *@ geno_param: reference to a struct that contains all information concerning the genotyping phase
 *@ kmer
 *@ qual_kmer
 *@ offset_value
 *@
 */

/*
 struct snp_aux_table
 {
 uint32_t pos_list[AUX_TABLE_COLS];
 snp_info snp_list[AUX_TABLE_COLS];
 }__attribute__((packed));
 */

/**
 * @geno_param: reference to a struct that contains all information concerning the genotyping phase
 * @kmer
 * @qual_kmer
 * @snpSD
 * @offset_value
 */

static void index_table_clear_index(IndexTable *index_table, uint32_t index)
{
	index_table->table[index % INDEX_TABLE_SLOT_COUNT].count = 0;
}

static void index_table_clear(IndexTable *index_table)
{
	index_table->best = NULL;
	index_table->ambiguous = false;

	for (size_t i = 0; i < INDEX_TABLE_SLOT_COUNT; i++)
	{
		index_table_clear_index(index_table, i);
	}
}

/*
 * REMEMBER:
 * 	-index: Position reads within the reference genome
 * 	-kmer_pos: k-mer position in the reference genome
 */
static void improved_index_table_add(IndexTable *index_table, uint32_t index, uint32_t kmer_pos,
										unordered_map<uint32_t, unordered_set<uint32_t>> & index_2_kmer_pos_set, bool is_neighbor)
{
	if (is_neighbor)
	{
		/*
		 * We need to make sure that at least one origin k-mer support current position.
		 * If previously I didn't find a k-mer that supported that position within the reference genome
		 * which is not a k-mer belonging to the neighborhood don't improved index table.
		 */
		if (index_2_kmer_pos_set.find(index) == index_2_kmer_pos_set.end())
		{
			if (verbose_map > 2)
				std::cerr << "(1339)exit" << std::endl;
			return;
		}
	}

	/*
	 * INDEX_TABLE_SLOT_COUNT is 1009
	 * Remember: index is the position of the read inside the reference genome
	 */
	size_t slot_index = index % INDEX_TABLE_SLOT_COUNT;

	/*
	 * struct IndexTableSlot has:
	 * 	-IndexTableEntry entries[INDEX_TABLE_ENTRY_DEPTH]
	 * 	-unsigned short count: number element inside entries
	 */
	IndexTableSlot *slot = &index_table->table[slot_index];

	/*
	 * struct IndexTableEntry has:
	 *  -uint32_t index
	 *  -uint8_t freq: number k-mer that support that position
	 */
	IndexTableEntry *target = NULL;

//Iteratively, consider all the elements of array entries
	for (unsigned short i = 0; i < slot->count; i++)
	{
		IndexTableEntry *e = &slot->entries[i];

//Check if an element with the same index is present inside array entries
//Ovviamente tutti gli elementi all'interno di entrie sono diversi,
//trovato uno posso fare break.
		if (e->index == index)
		{
			++e->freq;
			target = e;
			break;
		}
	}
	if (target == NULL)
	{
		IndexTableEntry i
		{index, 1};

		//Forse queste due istruzioni si potevano fare in un'unica operazione
		slot->entries[slot->count] = i;
		target = &slot->entries[slot->count];
		++slot->count;
	}

//Keep track that there is a k-mer that supports index position
	index_2_kmer_pos_set[index].insert(kmer_pos);

// if no two kmer pos, do not consider it as a match
//remember that the first can't be one that belongs to the neighborhood
	if (index_2_kmer_pos_set[index].size() <= 1)
		return;

	if (index_table->best == NULL)
	{
		index_table->best = target;
		index_table->ambiguous = false;
	}
	else if (target == index_table->best)
	{
		index_table->ambiguous = false;
	}
//I find another position within the reference genome having the same number of k-mer that support it.
	else if (target->freq == index_table->best->freq)
	{
		index_table->ambiguous = true;
	}
	/*
	 * If I found an index within the reference genome that is supported by
	 * more k-mer than the current best, I update best
	 */
	else if (target->freq > index_table->best->freq)
	{
		index_table->best = target;
		index_table->ambiguous = false;
	}

	if (verbose_map > 2)
	{
		std::cerr << "\n\n**Current-Best: " << unsigned(index_table->best->index) << std::endl;
		std::cerr << "Ambiguous: " << index_table->ambiguous << std::endl;
		std::cerr << "Frequence: " << unsigned(index_table->best->freq) << std::endl;
	}
}

void index_table_add(IndexTable *index_table, uint32_t index)
{
	size_t slot_index = index % INDEX_TABLE_SLOT_COUNT;
	IndexTableSlot *slot = &index_table->table[slot_index];
	IndexTableEntry *target = NULL;

	for (unsigned short i = 0; i < slot->count; i++)
	{
		IndexTableEntry *e = &slot->entries[i];
		if (e->index == index)
		{
			++e->freq;
			target = e;
			break;
		}
	}

	if (target == NULL)
	{
#if DEBUG
		assert(slot->count < INDEX_TABLE_ENTRY_DEPTH);
#endif
		IndexTableEntry i
		{index, 1};
		slot->entries[slot->count] = i;
		target = &slot->entries[slot->count];
		++slot->count;
	}

	if (index_table->best == NULL)
	{
		index_table->best = target;
		index_table->ambiguous = false;
	}
	else if (target == index_table->best)
	{
		index_table->ambiguous = false;
	}
	else if (target->freq == index_table->best->freq)
	{
		index_table->ambiguous = true;
	}
	else if (target->freq > index_table->best->freq)
	{
		index_table->best = target;
		index_table->ambiguous = false;
	}
}
