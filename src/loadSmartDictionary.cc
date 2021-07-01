#include <assert.h>
#include "utility.h"
#include "geno.h"
#include "loadSmartDictionary.h"
#include "pileup.h"
#include "lib_aln_inexact_matching.h"
#include <inttypes.h>
#include <unordered_map>

//only for debug
int compleat_dictionary_verbose = 0;
int load_memory_verbose = 0;

void checkValue(const uint64_t snp_dict_size, const uint64_t numberCompleateKmer);

bwaidx_t* loadIndex(std::string path)
{
	string file = get_file_from_path(path);
	vector < string > nameFileWhitExtension = split_vargeno(file, '.');
	string nameFile = "output/" + nameFileWhitExtension[0];
	return lib_aln_idx_load(nameFile.c_str());
}

//void PopulateAuxSnp(ifstream &snp_Smartdict, snp_aux_table* snp_aux_table, uint32_t index, uint64_t size_table, uint8_t num_of_pos);

/*
 * Load in RAM the dictionary. Remember that the dictionary consists of three data structures:
 * an array(for single k-mer), an unordered map(for r.c. of complete k-mer), and an array of *snp_aux_table
 * (k-mer that occur within the reference genome several times)
 *
 * @geno_param
 */
SNPSmartDictionary* load_smartSNP_Dict(param_creation_geno* geno_param)
{
	//Open a channel to read smart dictionary SNP
	ifstream snp_Smartdict_file(geno_param->snp_dict_filename, ios::in | ios::binary);
	assert(snp_Smartdict_file);

	//A complete k-mer is a k-mer that has both direct and r.c. inside reference genome.
	//(in any case in the re-assembly we have one of the two).
	//A single k-mer is a k-mer that inside the reference genome appears direct or r.c.

	// Header of dictionary is composed by this 5 field:
	// Size of the dictionary(64bit)
	//Size of auxiliary table(64bit)
	//Value of k(8 bit)
	//Lenth of reassembly(64 bit)
	//NUmber of k-mer compleate(64 bit)

	//Read the size of the dictionary of the SNPs and of the auxiliary table
	const size_t snp_dict_size = read_uint64_t(snp_Smartdict_file);
	const size_t snp_aux_table_size = read_uint64_t(snp_Smartdict_file);

	//Get the value of k
	const uint8_t k = read_uint8_t(snp_Smartdict_file);
	geno_param->k = k;

	//Get len of re-assembly
	const uint64_t len_reassembly = read_uint64_t(snp_Smartdict_file);	//len_reassembly = 653186798;

	//Get number of compleate k-mer inside snp smart dictionary
	const uint64_t numberCompleateKmer = read_uint64_t(snp_Smartdict_file);

	//checkValue(snp_dict_size, numberCompleateKmer);

	//PrincipalDictionaty* snp_smart_dictionary = principalDictionary_init(snp_dict_size);
	Entry_dictionary_snp * snp_smart_dictionary = (Entry_dictionary_snp*) calloc(len_reassembly, sizeof(Entry_dictionary_snp));
	assert(snp_smart_dictionary);

	// Pointer to the snp_aux_table vector that contains the start position(and other info such as where the SNP is
	// inside k-mer) for k-mer that occur within the reference genome several times.

	snp_aux_table* snp_auxiliary_table = new snp_aux_table[snp_aux_table_size];
	assert(snp_auxiliary_table);

	// If a k-mer is complete, r.c. is stored inside hash map. n_compleat_kmer must be a power of 2
	TotaleKmerSnp *reverse_total_kmer_snp = total_kmer_snp_struct_init(numberCompleateKmer);	//16777216	<-ultimo//4194304
	assert(reverse_total_kmer_snp);

	fprintf(stderr, "SIZE DIZIONARIO SNP: %zu \n", snp_dict_size);
	fprintf(stderr, "SIZE AUXILIARY TABLE: %zu \n", snp_aux_table_size);
	fprintf(stderr,"k : %" PRIu8 "\n", k);
	fprintf(stderr,"len_reassembly : %" PRIu64 "\n", len_reassembly);
	fprintf(stderr,"Number of compleate k-mer : %" PRIu64 "\n", numberCompleateKmer);

	//Read every record of snp smart dictionary
	for (size_t i = 0; i < snp_dict_size; i++)
	{
		//For each k-mer inside snp.genome, recover from dictionary the following information
		const uint32_t pos_kmer_reassembly = read_uint32_t(snp_Smartdict_file);
		const uint32_t pos_ref_genome = read_uint32_t(snp_Smartdict_file);	//Posizione all'interno del genoma di riferimento della k-mer "originale"
		const uint8_t flag = read_uint8_t(snp_Smartdict_file);
		const uint8_t snp = read_uint8_t(snp_Smartdict_file);			//pos all'interno della k-mer+riferimento
		const uint8_t alt_base = read_uint8_t(snp_Smartdict_file);
		const uint8_t ref_freq = read_uint8_t(snp_Smartdict_file);
		const uint8_t alt_freq = read_uint8_t(snp_Smartdict_file);

		/*
		 * snp is an 8-bit field, the 3 least significant bits define the base SNP, the most significant 5 bits
		 * the position of the SNP inside the k-mer
		 */
		const unsigned snp_info_ref = SNP_INFO_REF(snp);
		const uint8_t snp_info_pos = SNP_INFO_POS(snp);

		//Check if this k-mer is direct or r.c. inside the re-assembly compared to how it occurs in the reference
		uint8_t bit_rev = (flag & 0x01);

		uint8_t bit_compleat = ((flag & 0x02) >> 1);

		uint8_t bit_flag_ambig_for = ((flag & 0x4) >> 2);

		uint8_t bit_flag_ambig_rev = ((flag & 0x8) >> 3);

		uint8_t bit_flag_ambig = 0;

		if (bit_rev == 0)
			bit_flag_ambig = bit_flag_ambig_for;
		else
			bit_flag_ambig = bit_flag_ambig_rev;

		if ((bit_compleat == 0) || ((bit_compleat == 1) && (bit_rev == 0)))
		{
			snp_smart_dictionary[pos_kmer_reassembly].pos_snp_genome = pos_ref_genome;
			snp_smart_dictionary[pos_kmer_reassembly].is_reverse = bit_rev;
			snp_smart_dictionary[pos_kmer_reassembly].is_compleate = bit_compleat;
			snp_smart_dictionary[pos_kmer_reassembly].isDirectAmb = bit_flag_ambig_for;
			snp_smart_dictionary[pos_kmer_reassembly].isRevAmb = bit_flag_ambig_rev;
			snp_smart_dictionary[pos_kmer_reassembly].altSnpPosInsideKmer = snp_info_pos;
		}
		else
			//insert info inside unordered map
			total_kmer_snp_struct_add(reverse_total_kmer_snp, pos_kmer_reassembly, pos_ref_genome, snp_info_pos);

		//Add info pileup table

		// i.e. if the reference base is A, C, G or T
		//if ((snp_info_ref & BASE_N) == 0 && pos_ref_genome != POS_AMBIGUOUS && bit_flag_ambig == FLAG_UNAMBIGUOUS)
		if ((snp_info_ref & BASE_N) == 0 && pos_ref_genome != POS_AMBIGUOUS && bit_flag_ambig == FLAG_UNAMBIGUOUS)
		{
			//Position of the SNPs within the reference genome(index for pileup table)
			const uint32_t snp_pos = pos_ref_genome + snp_info_pos;

			ptable_add(&geno_param->ptable, snp_pos, snp_info_ref, alt_base, ref_freq, alt_freq);
		}
	}

	for (int i = 0; i < snp_aux_table_size; i++)
	{
		snp_auxiliary_table[i].pos_list = new uint32_t[AUX_TABLE_COLS];
		snp_auxiliary_table[i].snp_list = new uint8_t[AUX_TABLE_COLS];

		for (size_t j = 0; j < AUX_TABLE_COLS; j++)
		{
			snp_auxiliary_table[i].pos_list[j] = read_uint32_t(snp_Smartdict_file);
			//snp_auxiliary_table[i].snp_list[j] = read_uint8_t(snp_Smartdict_file);

			uint8_t temp = read_uint8_t(snp_Smartdict_file);
			snp_auxiliary_table[i].snp_list[j] = SNP_INFO_POS(temp);

			const uint8_t alt_base = read_uint8_t(snp_Smartdict_file);
			const uint8_t ref_freq = read_uint8_t(snp_Smartdict_file);
			const uint8_t alt_freq = read_uint8_t(snp_Smartdict_file);
		}
	}

	snp_Smartdict_file.close();
	return new SNPSmartDictionary
	{snp_smart_dictionary, len_reassembly, snp_auxiliary_table, snp_aux_table_size, reverse_total_kmer_snp};
}
/*
 * For each k-mer stored in the temporary file, search it inside the re-assembly in order to complete the
 * pos_reassembly field (of each record) in the smart dictionary. Moreover, i update the field "flag" of each record
 * if a k-mer is compleate by set the second bit less significant to 1.
 *
 * @path_smartDict: Dictionary path
 * @path_reassembly_snp: Re-assembly path
 * @snpTemp_filename: Path where k-mer is stored
 */
void complete_smart_dict(std::string path_smartDict, std::string path_reassembly_snp, std::string snpTemp_filename)
{
	/*
	 * Counter of the number of complete k-mer(k-mer present forward and r.c. inside reference genome modify)
	 * This information will store inside the header of dictionaty
	 */
	uint64_t cnt_compleat_kmer = 0;

	//Keeps track of the number of the current record analyzed inside dictionary
	uint64_t currentNumberRecord = 0;

	//Contains, iteratively, the contents of the flag field stored in the dictionary at the record @currentNumberRecord.
	uint8_t currentFlag = 0;

	//Counter of the number of complete k-mer
	uint64_t num_kmer_correct_compleate = 0;

	//Open a stream to read k-mer store inside snpTemp
	ifstream snp_kmer_temp(snpTemp_filename);

	std::string header;

	//K-mer inside tempFIle
	std::string kmer_from_temp_file;

	//Load in memory index
	bwaidx_t* index_reassembly = loadIndex(path_reassembly_snp.c_str());
	assert(index_reassembly);

	//Used to write the re-assembly position inside smart dictionary and to update field flag inside dictionary
	ofstream smartDictionary_write(path_smartDict, ios::binary | ios::out | ios::in);

	/*
	 * The first 16 bytes of the dictionary store the size of the dictionary and auxiliary table.
	 * In addition, another byte is used to store the value of k.
	 * Another 8 bytes are dedicated to storing the length of the re-assembly
	 * Finally, another 8 bytes are dedicated to storing the length of the re-assembly
	 * So, I move BYTE_HEADER_SMART_DICT = 33 (byte) to skip the header
	 *
	 */
	smartDictionary_write.seekp(BYTE_HEADER_SMART_DICT, ios_base::beg);

	//Used to read the flag field of each record inside the dictionary
	ifstream smartDictionary_read(path_smartDict, ios::binary);

	//Point to the flag field of the first record inside the dictionary(jump pos_reassembly and pos_genome/index aux table)
	smartDictionary_read.seekg(BYTE_HEADER_SMART_DICT + 2 * BYTE_STORE_POS, ios_base::beg);

	//uint32_t* checkCompleate = (uint32_t*) calloc(6000, sizeof(uint32_t));

	//When identify the relative position insidere re-assembly of a k-mer,
	unordered_map < uint32_t, uint64_t > checkCompleate;

	while (std::getline(snp_kmer_temp, header))
	{
		std::getline(snp_kmer_temp, kmer_from_temp_file);

		currentNumberRecord++;

		//Set the local variable flag with the value of the field flag inside the dictionary
		smartDictionary_read.read(reinterpret_cast<char *>(&currentFlag), sizeof(uint8_t));

		if (compleat_dictionary_verbose > 2)
			std::cerr << "\nThe current k-mer to be found in re-assembly is: " << kmer_from_temp_file << std::endl;

		//Number of hit found by lib_aln_bound_backtracking(always 1) mm = 0 by default
		hit_reassembly* result = lib_aln_bound_backtracking_record(index_reassembly, kmer_from_temp_file.c_str());
		/*
		 * Initially "currentFlag" is equal to 1 if current k-mer is ambiguous, while is equal to 0 if the current k-mer is no ambiguous.
		 * Simply copy this value inside variable currentFlag_bitAmbig; this variable is useful is k-mer is compleate
		 */
		uint8_t currentFlag_bitAmbig = currentFlag;

		//La k-mer associata al record corrente nel dizionario è presente reverse nel reassembly
		if (result->isRev == true)
		{
			//Nel caso in cui è reverse, devo spostare di 4 bit a sinistra
			currentFlag = currentFlag << 3;
			/*
			 * Set the first bit of flag to 1 => this k-mer is present r.c. in the reassembly
			 * respect to the reference genome
			 */
			currentFlag = (currentFlag | 0x01);
		}
		else
		{
			//Semplicemente sposto il bit ambiguità alla posizione 3 in quanto tale k-mer è stata trovata diretta
			currentFlag = currentFlag << 2;
		}

		//Position of k-mer inside reassembly
		uint32_t posf = result->positions_to_ref[0];

		if (compleat_dictionary_verbose > 2)
			std::cerr << "The position in the reassembly results: " << posf << std::endl;

		//Check if previusly i found this k-mer directo o r.c
		if (checkCompleate.find(posf) == checkCompleate.end())
			checkCompleate[posf] = currentNumberRecord;
		else
		{
			//Get the previous number of the other k-mer that is
			uint32_t n_rec = checkCompleate.at(posf);

			//This means that the k-mer is compleate
			cnt_compleat_kmer++;

			/**
			 * N.B- If the current k-mer is forward, means that the first one war r.c.
			 * Otherwise, if this k-mer is r.c, this means that this first one was forward.
			 */

			//Tengo traccia che la k-mer è completa, il secondo bit meno significativo del campo flag corrente viene posto a 1
			currentFlag = (currentFlag | 0x2);

			//Recupero il campo flag della k-mer vista prima
			uint8_t previousFlag_bitAmbig = 0;

			//Recupero il campo flag della k-mer vista prima
			uint8_t previousFlag = 0;

			//Used to write the re-assembly position inside smart dictionary and to update field flag inside dictionary
			ofstream smartDictionary_write_compleate(path_smartDict, ios::binary | ios::out | ios::in);

			//Used to read the flag field of each record inside the dictionary
			ifstream smartDictionary_read_compleate(path_smartDict, ios::binary);

			smartDictionary_read_compleate.seekg(BYTE_HEADER_SMART_DICT + (n_rec - 1) * BYTE_REC_SNP + 8, ios_base::beg);

			//Set the local variable flag with the value of the field flag inside the dictionary
			smartDictionary_read_compleate.read(reinterpret_cast<char *>(&previousFlag), sizeof(uint8_t));

			if (result->isRev == true)
			{
				//Il flag ambiguità della precedente è memorizzato in pos 0x4 = 00000100
				previousFlag_bitAmbig = (previousFlag & 0x4) >> 2;
				//Se questa è reverse aggiungo questa info alla precedente che era diretta
				previousFlag = (previousFlag | (currentFlag_bitAmbig << 3));
				currentFlag = currentFlag | (previousFlag_bitAmbig << 2);
			}
			else
			{
				// La k-mer incontrata prima è la r.c, conserva il bit amb alla pos 4 = 0x8 = 00001000
				previousFlag_bitAmbig = (previousFlag & 0x8) >> 3;

				/*
				 * Se invece questa è quella diretta, precedentemente ho visto l'inversa e quindi
				 * ad essa aggiungo l'informazione della diretta
				 */
				previousFlag = (previousFlag | (currentFlag_bitAmbig << 2));

				//Qella di prima è r.c.
				currentFlag = currentFlag | (previousFlag_bitAmbig << 3);
			}

			if (posf == 320587696)
			{
				uint8_t bit_rev = (previousFlag & 0x01);

				uint8_t bit_compleat = ((previousFlag & 0x02) >> 1);

				uint8_t bit_flag_ambig_for = ((previousFlag & 0x4) >> 2);

				uint8_t bit_flag_ambig_rev = ((previousFlag & 0x8) >> 3);

			fprintf(stderr,"\nPREVIEUS FLAG-ORIGINAL: %"PRIu32 "\n", posf);
			fprintf(stderr,"pos_kmer_reassembly: %"PRIu32 "\n", posf);
			fprintf(stderr,"Bit rev : %" PRIu8 "\n", bit_rev);
			fprintf(stderr,"Bit bit_compleat : %" PRIu8 "\n", bit_compleat);
			fprintf(stderr, "Ambig For : %" PRIu8 "\n", bit_flag_ambig_for);
			fprintf(stderr, "Ambig rev : %" PRIu8 "\n", bit_flag_ambig_rev);
		}
		//The k-mer is compleate
		previousFlag = (previousFlag | 0x2);

		smartDictionary_write_compleate.seekp(BYTE_HEADER_SMART_DICT + (n_rec - 1) * BYTE_REC_SNP + 8, ios_base::beg);
		smartDictionary_write_compleate.write((char *) (&previousFlag), sizeof(uint8_t));

		smartDictionary_write_compleate.flush();
		smartDictionary_write_compleate.close();
		smartDictionary_read_compleate.close();

		if (posf == 320587696)
		{
			uint8_t bit_rev = (previousFlag & 0x01);

			uint8_t bit_compleat = ((previousFlag & 0x02) >> 1);

			uint8_t bit_flag_ambig_for = ((previousFlag & 0x4) >> 2);

			uint8_t bit_flag_ambig_rev = ((previousFlag & 0x8) >> 3);

		fprintf(stderr,"\nPREVIEUS FLAG-AFTER: %"PRIu32 "\n", posf);
		fprintf(stderr,"pos_kmer_reassembly: %"PRIu32 "\n", posf);
		fprintf(stderr,"Bit rev : %" PRIu8 "\n", bit_rev);
		fprintf(stderr,"Bit bit_compleat : %" PRIu8 "\n", bit_compleat);
		fprintf(stderr, "Ambig For : %" PRIu8 "\n", bit_flag_ambig_for);
		fprintf(stderr, "Ambig rev : %" PRIu8 "\n", bit_flag_ambig_rev);
	}

}

//Write the position inside the re-assembly in the dictionary
smartDictionary_write.write((char *) (&(posf)), sizeof(uint32_t));

//I move to write in the flag field
smartDictionary_write.seekp(BYTE_STORE_POS, ios_base::cur);

//Write the field flag
smartDictionary_write.write((char *) (&currentFlag), sizeof(uint8_t));

//Set the pointer for write the next field pos reassembly.
smartDictionary_write.seekp(4, ios_base::cur);
//Set the pointer for read the next "field flag"
smartDictionary_read.seekg(COMPLETE_SNP_DICT, ios_base::cur);

if (posf == 320587696)
{
	uint8_t bit_rev = (currentFlag & 0x01);

	uint8_t bit_compleat = ((currentFlag & 0x02) >> 1);

	uint8_t bit_flag_ambig_for = ((currentFlag & 0x4) >> 2);

	uint8_t bit_flag_ambig_rev = ((currentFlag & 0x8) >> 3);

fprintf(stderr,"pos_kmer_reassembly: %"PRIu32 "\n", posf);
fprintf(stderr,"Bit rev : %" PRIu8 "\n", bit_rev);
fprintf(stderr,"Bit bit_compleat : %" PRIu8 "\n", bit_compleat);
fprintf(stderr, "Ambig For : %" PRIu8 "\n", bit_flag_ambig_for);
fprintf(stderr, "Ambig rev : %" PRIu8 "\n", bit_flag_ambig_rev);
}

//free the memory from result
lib_aln_hit_reassembly_destroy(result);

num_kmer_correct_compleate++;

if ((num_kmer_correct_compleate % 100000000) == 0)
{
std::cerr << "Num read correct compleate is 100000000" << debugProgress() << std::endl;
}
}		//End map every k-mer inside temp file

//std::cerr << "Num read correct compleate is: " << (unsigned) num_kmer_correct_compleate << std::endl;
std::cerr << "Numner of  compleate k-mer is: " << (unsigned) cnt_compleat_kmer << std::endl;

//free (checkCompleate);
smartDictionary_write.seekp(0);
smartDictionary_write.seekp(8 + 8 + 1 + 8, ios_base::beg);
smartDictionary_write.write((char *) (&cnt_compleat_kmer), sizeof(uint64_t));
//Free memory for index
lib_aln_idx_destroy(index_reassembly);

smartDictionary_write.flush();
smartDictionary_write.close();
smartDictionary_read.close();
}

void checkValue(const uint64_t snp_dict_size, const uint64_t numberCompleateKmer)
{
if (snp_dict_size > POW_2_32)
{
fprintf(stderr, "Error: SNP dictionary is too large (limit: %lu 32-mers)\n", POW_2_32);
exit (EXIT_FAILURE);
}

if (snp_dict_size <= 0)
{
std::cerr << "Error: number k-mer inside dictionary can't be < 0" << std::endl;
exit (EXIT_FAILURE);
}

if (numberCompleateKmer < 0)
{
//Thin number can be equal to 0
fprintf(stderr, "Error: compleate k-mer inside dictionary can't be < 0");
exit (EXIT_FAILURE);
}

}
