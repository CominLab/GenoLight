#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <inttypes.h>
#include <algorithm>
#include "dictgen.h"
#include "utility.h"
#include "ParseFasta.hpp"
#include <ctime>
#include <assert.h>

//Only for debug
uint8_t verbose_dictgen = 0;

std::string decode_kmer_vargeno(uint64_t e, uint8_t k_length);
bool ExtractFrequency(std::string info, float& f1, float& f2);

/**
 * Extract the input parameters present in the command entered by the user
 *
 * @argc
 * @argv
 */
param_creation_dict vargeno_init(int argc, char* argv[])
{
	param_creation_dict dict_param;
	dict_param.k = 32;

	std::string prefix;

	//read the parameters
	for (int i = 2; i < argc; ++i)
	{
		if (strcmp(argv[i], "-k") == 0)
		{
			//k-mer length
			++i;
			dict_param.k = strtoull(argv[i], nullptr, 10);
		}
		if (strcmp(argv[i], "-s") == 0)
		{
			++i;
			dict_param.snp_filename = argv[i];
		}
		else if (strcmp(argv[i], "-r") == 0)
		{

			++i;
			dict_param.ref_filename = argv[i];
		}
		else if (strcmp(argv[i], "-p") == 0)
		{

			++i;
			prefix = argv[i];
		}
		else
		{
			void print_help_building_dictionaries();
			i = argc;
		}
	}

	create_default_out_folder("output/");

	if (prefix.size() == 0)
		prefix = get_name_from_genome_path(dict_param.ref_filename);

	dict_param.chrlens_filename = "output/" + prefix + ".chrlens";

	dict_param.path_info_file = "output/" + prefix + ".info";
	dict_param.path_ref_genome_index = "output/" + prefix;

	std::string s = get_file_from_path(dict_param.snp_filename);
	std::vector < std::string > col = split_vargeno(s, '.');

	//Devo dare una estensione in modo che sia semmetrico all'estensione .fa quando faccio l'indicizzazione di entrambi
	dict_param.path_reassembly_snp = "output/" + prefix + "Reassembly" + col[col.size() - 2] + ".txt";
	dict_param.snp_dict_filename = "output/" + prefix + "SmartDict" + col[col.size() - 2];
	dict_param.snpTemp_filename = "output/" + prefix + "Temp" + ".txt";
	return dict_param;
}
/**
 * @e
 */
void decode_8(uint8_t e)
{
	int i;
	//printf("\nLa cod è: %" PRIu64 "\n", e);

	for (i = 0; i < 8; ++i)
	{
		switch (e & 0x03)
		{
			case 0:
			{
				printf("00");
			}
				break;
			case 1:
			{
				printf("01");
			}

				break;

			case 2:
			{
				printf("10");
			}
				break;
			case 3:
			{
				printf("11");
			}
				break;
		}
		e >>= 2;
	}
}

/*
 * For each k-mer inside vector kmers, it stores in temp file its literal encoding and in the dictionary a records
 * with all its information.
 *
 * @kmers: This vector contains all the k-mer extracted from the modified reference genome.
 * @kmers_len
 * @snp_dict_filename: Path to store the dictionary
 * @snpTemp_filename: Path to store k-mer snpTemp_filename
 * @k_length:
 */
void write_snp_kmers(struct snp_kmer_info *kmers, const size_t kmers_len, std::string snp_dict_filename, std::string snpTemp_filename,
						uint8_t k_length)
{
	//Open the flow to create the smartDictionary of the SNPs
	//std::ofstream snpdict_file(snp_dict_filename, std::ios::binary);

	FILE *snpdict_file = fopen(snp_dict_filename.c_str(), "wb");
	assert(snpdict_file);

	//Open the flow to create the temporary file whit the kmer containing a SNP

	std::ofstream snpTemp_file(snpTemp_filename);

	if (!snpTemp_file.good())
	{
		std::cerr << "(146-dictgen.cc)Error create: " << snpTemp_filename << std::endl;
		exit (EXIT_FAILURE);
	}

	/*
	 * NOTE:
	 * -In the dictionary binary file, first I store the dictionary and then the auxiliary table.
	 */
	std::vector < snp_aux_table_dictgen > aux_table(SNP_AUX_TABLE_INIT_SIZE);

	//Number of elements present in the auxiliary table
	uint64_t aux_table_count = 0;

	//I leave the space to store the size(rows) of the dictionary and of the auxiliary table
	const uint64_t placeholder = 0UL;

	serialize_uint64(snpdict_file, placeholder); /* dict size (rows) placeholder */
	serialize_uint64(snpdict_file, placeholder); /* aux table size (rows) placeholder */

	//write_uint64_t(snpdict_file, placeholder);
	//write_uint64_t(snpdict_file, placeholder);

	//I memorize the value of k
	serialize_uint8(snpdict_file, k_length);

	//I memorize the len of reassembly of snp
	serialize_uint64(snpdict_file, placeholder);

	//I memorize the number of compleate kmer
	serialize_uint64(snpdict_file, placeholder);

	//Number of k-mer stored in the dictionary
	uint64_t kmers_written = 0UL;

	//Indicates the index of the current k-mer that is stored
	size_t i = 0;

	/**
	 * Keep track of a few statistics:
	 *
	 * total_kmers: Number of k-mer extracted from reference genome modify whit snp
	 *
	 * unambig_kmers: Number of k-mer that have a number of occurrences of 1
	 *
	 * ambig_unique_kmers: Number of distinct k-mers having a multiplicity greater than 1
	 *
	 * ambig_total_kmers: Total number of k-mer ambiguous (ie that in place of the position they have in the pos field the index to access the auxiliary table)
	 *
	 * total_pos_ambig_kmers: k-mer having a greater multiplicity of AUX_TABLE_COLS
	 */

	const size_t total_kmers = kmers_len;
	size_t unambig_kmers = 0;
	size_t ambig_unique_kmers = 0;
	size_t ambig_total_kmers = 0;
	size_t total_pos_ambig_kmers = 0;

	std::cerr << (unsigned) kmers_len << std::endl;

	while (i < kmers_len)
	{
		//Current k-mer
		const kmer_t kmer = kmers[i].kmer;

		//I write the k-mer in the temporary text file
		//Convert k-mer to binary format into string format
		std::string kmer_decod = decode_kmer_vargeno(kmer, k_length);

		std::string record = ">Number k-mer: " + std::to_string(i) + "\n" + kmer_decod + "\n";
		snpTemp_file.write(record.c_str(), sizeof(char) * record.size());

		//Position of k-mer inside the reference genome
		const uint32_t pos = kmers[i].pos;

		/*
		 * snp_info is an 8-bit integer
		 * 	-5 bits: position of the SNP within the k-mer
		 * 	-3 bits: reference base(ie base present in the reference genome)
		 */
		const snp_info snp = kmers[i].snp;

		//base with which it is replaced reference base
		const uint8_t alt_base = kmers[i].alt_base;

		const uint8_t ref_freq = kmers[i].ref_freq;
		const uint8_t alt_freq = kmers[i].alt_freq;

		++i;

		//Check if this k-mer has more occurrence within the reference genome
		if (i < kmers_len && kmer == kmers[i].kmer)
		{
			++ambig_unique_kmers;
			++ambig_total_kmers;

			//Memorizzo subito la k-mer attuale
			aux_table[aux_table_count].pos_list[0] = pos;
			aux_table[aux_table_count].snp_list[0] = snp;
			aux_table[aux_table_count].alt_base[0] = alt_base;
			aux_table[aux_table_count].ref_freqs[0] = ref_freq;
			aux_table[aux_table_count].alt_freqs[0] = alt_freq;

			size_t k = 1;
			bool too_many_positions = false;

			/*
			 * I store(momentarily in RAM) all the others up to a maximum of 8
			 * As mentioned at the beginning, for each one I have to memorize all the information as:
			 * initial position, where is SNPs, what is the reference base etc.
			 */
			do
			{
				//AUX_TABLE_COLS is equal to 10
				if (k < AUX_TABLE_COLS)
				{
					aux_table[aux_table_count].pos_list[k] = kmers[i].pos;
					aux_table[aux_table_count].snp_list[k] = kmers[i].snp;
					aux_table[aux_table_count].alt_base[k] = kmers[i].alt_base;
					aux_table[aux_table_count].ref_freqs[k] = kmers[i].ref_freq;
					aux_table[aux_table_count].alt_freqs[k] = kmers[i].alt_freq;
					++k;
				}
				else
				{
					total_pos_ambig_kmers++;
					too_many_positions = true;
				}

				++ambig_total_kmers;
				++i;
			}
			while (i < kmers_len && kmer == kmers[i].kmer);

			/*
			 * If the current k-mer occurs more than 10 times within the reference genome,
			 * I store ONLY in the main dictionary the information that that k-mer will
			 * not be managed by setting POS to POS_AMBIGUOUS.
			 */
			if (too_many_positions)
			{
				serialize_uint32(snpdict_file, 0);
				serialize_uint32(snpdict_file, POS_AMBIGUOUS);
				serialize_uint8(snpdict_file, FLAG_AMBIGUOUS);
				serialize_uint8(snpdict_file, 0);  // SNP info
				serialize_uint8(snpdict_file, 0);  // SNP alt
				serialize_uint8(snpdict_file, 0);  // ref freq
				serialize_uint8(snpdict_file, 0);  // alt freq
			}
			else
			{
				//uint8_t flag = (k << 2);

				// fill in remainder with 0s
				while (k < AUX_TABLE_COLS)
				{
					aux_table[aux_table_count].pos_list[k] = 0;
					aux_table[aux_table_count].snp_list[k] = 0;
					aux_table[aux_table_count].alt_base[k] = 0;
					aux_table[aux_table_count].ref_freqs[k] = 0;
					aux_table[aux_table_count].alt_freqs[k] = 0;
					++k;
				}

				serialize_uint32(snpdict_file, 0);
				serialize_uint32(snpdict_file, aux_table_count);
				serialize_uint8(snpdict_file, FLAG_AMBIGUOUS); //Nel caso in cui ho un numero di duplicati inferiore a 10 memorizzo in flag il numero di tali occorrenze
				serialize_uint8(snpdict_file, 0);  // SNP info
				serialize_uint8(snpdict_file, 0);  // SNP alt
				serialize_uint8(snpdict_file, 0);  // ref freq
				serialize_uint8(snpdict_file, 0);  // alt freq
				++aux_table_count;
			}
		}
		else
		{
			//Case in which the k-mer is unique
			++unambig_kmers;

			serialize_uint32(snpdict_file, 0);
			serialize_uint32(snpdict_file, pos);
			serialize_uint8(snpdict_file, FLAG_UNAMBIGUOUS);
			serialize_uint8(snpdict_file, snp);
			serialize_uint8(snpdict_file, alt_base);
			serialize_uint8(snpdict_file, ref_freq);  // ref freq
			serialize_uint8(snpdict_file, alt_freq);  // alt freq
		}

		++kmers_written;
	}

	std::cerr << "\nStart write kmer aux table" << std::endl;

	//Finally, I store the auxiliary table in the smart dictionary file
	for (i = 0; i < aux_table_count; i++)
	{
		for (size_t j = 0; j < AUX_TABLE_COLS; j++)
		{
			serialize_uint32(snpdict_file, aux_table[i].pos_list[j]);
			serialize_uint8(snpdict_file, aux_table[i].snp_list[j]);
			serialize_uint8(snpdict_file, aux_table[i].alt_base[j]);
			serialize_uint8(snpdict_file, aux_table[i].ref_freqs[j]);
			serialize_uint8(snpdict_file, aux_table[i].alt_freqs[j]);
		}
	}

	rewind(snpdict_file);

	//Store the size of the dictionary and auxiliary table at the beginning of the file
	serialize_uint64(snpdict_file, kmers_written);
	serialize_uint64(snpdict_file, aux_table_count);

	fclose(snpdict_file);

	snpTemp_file.flush();
	snpTemp_file.close();

	std::cerr << "\nSNP Dictionary summary: " << std::endl;

	std::cerr << "Total k-mers: " << total_kmers << std::endl;
	std::cerr << "Unambig k-mers: " << unambig_kmers << std::endl;
	std::cerr << "Ambig unique k-mers: " << ambig_unique_kmers << std::endl;
	std::cerr << "Ambig total k-mers: " << ambig_total_kmers << std::endl;
	std::cerr << "Ambig pos ambig k-mer: " << total_pos_ambig_kmers << std::endl;

	std::cout << " End of dictionary creation SNPs. The number of record inside dictionary is: " << kmers_written << '\n';
}

/**
 * @p1
 * @p2
 */
int kmer_cmp(const void* p1, const void* p2)
{
	const kmer_t kmer1 = ((struct kmer_info *) p1)->kmer;
	const kmer_t kmer2 = ((struct kmer_info *) p2)->kmer;
	return (kmer1 > kmer2) - (kmer1 < kmer2);
}

/**
 * @c
 */
static char rev(const char c)
{
	switch (c)
	{
		case 'A':
		case 'a':
			return 'T';
		case 'C':
		case 'c':
			return 'G';
		case 'G':
		case 'g':
			return 'C';
		case 'T':
		case 't':
			return 'A';
		default:
			return 'N';
	}
}

/**
 * @base
 */
static uint8_t convertBaseBinary(char base)
{
	switch (base)
	{
		case 'A':
		case 'a':
			return 0x0;
		case 'C':
		case 'c':
			return 0x1;
		case 'G':
		case 'g':
			return 0x2;
		case 'T':
		case 't':
			return 0x3;
		default:
			return 'N';
	}
}

int snp_kmer_cmp(const void *p1, const void *p2)
{
	const kmer_t kmer1 = ((struct snp_kmer_info *) p1)->kmer;
	const kmer_t kmer2 = ((struct snp_kmer_info *) p2)->kmer;
	return (kmer1 > kmer2) - (kmer1 < kmer2);
}

static void sort_snp_kmers(struct snp_kmer_info *kmers, const size_t kmers_len)
{
	qsort(kmers, kmers_len, sizeof(*kmers), snp_kmer_cmp);
}

/*
 * `snp_file` in this case is a file in the UCSC SNP-txt format, consisting
 * of common SNP to be included in the dictionary.
 *
 * http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=varRep&hgta_track=snp141Common&hgta_table=snp141Common&hgta_doSchema=describe+table+schema
 * k=snp141Common&hgta_table=snp141Common&hgta_doSchema=describe+table+schema
 *
 * Be sure to filter any SNPs with abnormal conditions (e.g. inconsistent alleles).
 *
 * @genome
 * @snp_dict_filename
 * @snp_filename
 * @snpTemp_filename
 * @snp_locations
 *
 * @return
 */
void make_snp_dict_from_vcf(ParseFasta genome, std::string snp_dict_filename, std::string snp_filename, std::string snpTemp_filename,
							bool **snp_locations, size_t snp_locs_size)
{

#define CHROM_FIELD   0
#define INDEX_FIELD   1
#define REF1_FIELD    3
#define REF2_FIELD    3
#define ALT_FIELD     4
#define INFO_FIELD	  7
#define FREQS_FIELD   25

	*snp_locations = (bool*) malloc(snp_locs_size * sizeof(bool));
	assert(*snp_locations);
	memset(*snp_locations, false, snp_locs_size);

	std::size_t numberValidSNP = 0;

	std::cerr << debugProgress() << std::endl;

	std::size_t max_pos_snp = 0;

	//Read the list of SNPs
	std::ifstream snp_file(snp_filename);

	if (!snp_file.good())
	{
		std::cerr << "(466-dictgen.cc)Error opening: " << snp_filename << std::endl;
		exit (EXIT_FAILURE);
	}

	//Store name of chromosome
	std::string chrom_name;

	//Store a record of .vcf file contain information snp
	std::string line_string;

	//Split line_string
	std::vector < std::string > line_split;

	//counter of the number of k-mer extracted
	size_t kmers_len = 0;

	unsigned int start_index_chromosome = 1;  // 1-based

	//ho imposto io di avercelo
	/*
	 bool ref_has_chr = true;
	 if (genome[0].name[0] != 'c')
	 ref_has_chr = false;
	 */

	FastaRecord* chrom = nullptr;

	//Defines if the SNP's record contains the frequencies associated to the reference base and the alternative base.
	bool has_freq = true;

	size_t n_vnc = 0;
	while (std::getline(snp_file, line_string))
	{
		/*
		 * Discard all the records related to Meta-information
		 * (for example fileformat,fileDate, reference ..)
		 */
		if (line_string.empty())
			continue;
		if (line_string[0] == '#')
			continue;

		++n_vnc;
	}
	snp_file.clear();
	snp_file.seekg(0, std::ios::beg);

	size_t max_kmers_len = 32 * n_vnc;
	struct snp_kmer_info *snp_kmers = (struct snp_kmer_info*) malloc(max_kmers_len * sizeof(*snp_kmers));
	assert(snp_kmers);

	//Initial and final instant of an operation that wants to calculate the duration

	clock_t begin, end;
	double time_spent;
	begin = clock();

	while (std::getline(snp_file, line_string))
	{
		/*
		 * Discard all the records related to Meta-information
		 * (for example fileformat,fileDate, reference ..)
		 */
		if (line_string.empty())
			continue;

		if (line_string[0] == '#')
			continue;

		std::vector < std::string > columns = split_vargeno(line_string, '\t');

		chrom_name = columns[CHROM_FIELD];

		if (columns[CHROM_FIELD][0] != 'c')
		{
			chrom_name = "chr" + columns[CHROM_FIELD];
		}

		// copy chromosome name into an independent buffer
		size_t chrom_index = stoi(columns[INDEX_FIELD]) - 1; // 1-based to 0-based coordinate

		//Read REF base
		const char ref_base = toupper(columns[REF1_FIELD][0]);

		//Encode REF base
		const unsigned ref_base_u = encode_base(ref_base);

		//BASE_X mean any other base other than A, C, G, T, N
		if (ref_base_u == BASE_X)
			continue;

		/* check if reference sequences are 1 base long */

		if (columns[REF1_FIELD].size() != 1)
			continue;

		if (columns[ALT_FIELD].size() != 1)
			continue;

		/*
		 * Questa verifica non è stupida, nel file .vcf possiamo avere più record consecutivi
		 * relativi a SNPs dello stesso cromosoma, in questo caso non ha senso che lo vada a ricercare
		 * tutte le volte.
		 * chrom == NULL per la prima volta
		 * chrom_name.compare(chrom->name) != 0 per quelle successiva
		 */

		/*
		 * In the .vcf file we can have multiple consecutive records related to SNPs of the same chromosome,
		 * for this reason it makes no sense to look for the chromosome several times.
		 */
		if ((numberValidSNP == 0) || chrom_name.compare(chrom->name) != 0)
		{
			chrom = genome.find_seq_by_name(chrom_name, start_index_chromosome);

			if (chrom == nullptr)
			{
				std::cerr << "\n[Error] chromosome name " << chrom_name
				<< " in VCF file not found in reference genome FASTA file.\n Usually this is because the FASTA file "
				"has chromosome name as \"chr1\" while the VCF file has chromosome name as \"1\" without the \"chr\"" << std::endl;
				continue;
			}
		}

		uint32_t seqChromSize = chrom->sequence.size();

		//Position within the CHROMOSOME where the SNP is located
		const unsigned int index = std::stoi(columns[INDEX_FIELD]) - 1;  // Convert from 1-based to 0-based

		//I check that the reference base present in the .vcf file corresponds to the base present in index position within the reference genome
		if (index >= seqChromSize)
		{
			std::cerr << "Error: The length of the chromosome " << chrom->name << " is minor than the pos of current snp " << index << std::endl;
			exit (EXIT_FAILURE);
			continue;
		}
		if (toupper(chrom->sequence[index]) != ref_base)
		{
			std::cerr << "Error: Mismatch found between reference sequence and SNP file at 0-based index " << index << " in " << chrom->name
			<< std::endl;
			std::cerr << "Reference base is: " << ref_base << std::endl;
			std::cerr << "genome is: " << chrom->sequence[index] << std::endl;
			exit (EXIT_FAILURE);
			continue;
		}

		/*
		 * I check if the index is too small or too large
		 * for every snp I have to make sure I can extract 32 k-mer)
		 */
		if (index < 32 || (index + 32) > seqChromSize)
			continue;

		/*
		 * We should only process bi-allelic SNPs.
		 * A biallelic site is a specific locus in a genome that contains two observed alleles, counting the reference as one,
		 * and therefore allowing for one variant allele.
		 */
		const bool neg = false;

		const char a1 = ref_base;
		const char a2 = toupper(columns[ALT_FIELD][0]);

		if (!(a1 == 'A' || a1 == 'C' || a1 == 'G' || a1 == 'T'))
			continue;
		if (!(a2 == 'A' || a2 == 'C' || a2 == 'G' || a2 == 'T'))
			continue;

		if (a1 != ref_base && a2 != ref_base)
		{
			//fprintf(stderr, "no ref base\n");
			continue;
		}

		//Is the position of the SNP inside the REFERENCE GENOME
		const size_t loc = start_index_chromosome + index;

		if (loc >= snp_locs_size)
		{
			*snp_locations = (bool*) realloc(snp_locations, (loc + 1) * sizeof(bool));
			assert(*snp_locations);
			memset(*snp_locations + snp_locs_size, false, loc - snp_locs_size + 1);
			snp_locs_size = loc + 1;
		}
		(*snp_locations)[loc] = true;

		if (max_pos_snp < loc)
			max_pos_snp = loc;

		std::string info = columns[INFO_FIELD];

		/*
		 * I assign default values to the frequency of the reference allele and to the alternative allele.
		 * This is because it is said that in the .vcf file this information is reported.
		 * If they are present now I'm going to recover them.
		 */
		float freq1 = 0.5;
		float freq2 = 0.5;

		uint8_t freq1_enc = 0;
		uint8_t freq2_enc = 0;

		if (ExtractFrequency(info, freq1, freq2))
		{
			freq1_enc = (uint8_t)(freq1 * 0xff);
			freq2_enc = (uint8_t)(freq2 * 0xff);
		}

		char *alt_p = const_cast<char*>(columns[ALT_FIELD].c_str());

		const char alt = neg ? rev(toupper(*alt_p)) : toupper(*alt_p);

		if (alt == ref_base || !(alt == 'A' || alt == 'C' || alt == 'G' || alt == 'T'))
			continue;

		//const std::string seq = chrom->sequence;

		//If the kmer before pos has N, then skip current SNP
		bool kmer_had_n;

		//It's the last k-mer that doesn't include the SNP
		kmer_t kmer = encode_kmer(&chrom->sequence[index - 32], kmer_had_n);

		//std::cerr << "LAST: " << decode_kmer_vargeno(kmer, dict_param->k) << std::endl;
		//std::cerr << "ALT: " << alt << std::endl;

		if (kmer_had_n)
			goto end;

		for (int i = 0; i < 32; i++)
		{
			if ((chrom->sequence[index + i] == 'N') || (chrom->sequence[index + i] == 'n'))
				goto end;
		}

		//std::cerr << "kmer"<<(unsigned)kmer << std::endl;

		for (unsigned int i = 0; i < 32; i++)
		{
			/*
			 * When i = 0, the condition results to be false,
			 * so next_base is equal to alt.
			 *
			 * Next the condition is always true so next_base in
			 * equal to chrom->sequence[index + i]
			 */
			const char next_base = (i ? chrom->sequence[index + i] : alt);
			//std::cerr << "NEXT BASE: " << next_base << std::endl;

			/*
			 * Once that the alternative base is "inside" k-mer, making all
			 * the possible rotations in k-mer we obtain the k-mer desired.
			 */
			kmer = shift_kmer(kmer, next_base);
			//std::cerr << "SHIFT K-MER: " << decode_kmer_vargeno(kmer, dict_param->k) << std::endl;

			//Remember that the kmer variable contains the encoding of the current k-mer.
			//Questa mi serve per memorizzare tale k-mer sul dizionario
			snp_kmers[kmers_len].kmer = kmer;

			//std::string temp = decode_kmer_vargeno(kmer, dict_param->k);
			//std::cerr << temp << std::endl;

			//Index is the position inside the genome of the SNPs
			snp_kmers[kmers_len].pos = start_index_chromosome + index - 32 + 1 + i;

			//Temporarily, I also memorize the alternative basis for simplicity.
			//Before to store, convert in bunary format{00,01,10,11}
			snp_kmers[kmers_len].alt_base = convertBaseBinary(alt);

			/*
			 * The snp field is of type snp_info (a 8-bit unsigned integer).
			 * 	-5 bits encode the position of the SNPs within the k-mer
			 * 	-3 bits for the reference base
			 */
			snp_kmers[kmers_len].snp = SNP_INFO_MAKE(32 - 1 - i, ref_base_u);

			//std::cerr << "SHIFT K-MER: " << decode_kmer_vargeno(kmer, 32) << std::endl;
			//printf("kmer %" PRIu64 "\n", kmer);
			//printf("pos: %" PRIu32 "\n", snp_kmers[i].pos);

			//I memorize the two frequencies
			snp_kmers[kmers_len].ref_freq = freq1_enc;
			snp_kmers[kmers_len].alt_freq = freq2_enc;
			kmers_len++;
		}

		++numberValidSNP;
		end: continue;
	}
#undef CHROM_FIELD
#undef INDEX_FIELD
#undef REF1_FIELD
#undef REF2_FIELD
#undef ALT_FIELD
#undef INFO_FIELD
#undef FREQS_FIELD


	std::cerr << "Number SNP:" << (unsigned)numberValidSNP << std::endl;
	snp_kmers = (struct snp_kmer_info*) realloc(snp_kmers, kmers_len * sizeof(*snp_kmers));

	sort_snp_kmers(snp_kmers, kmers_len);

	//Create the binary dictionary and the temporary text file where the literals are stored in the literal form
	write_snp_kmers(snp_kmers, kmers_len, snp_dict_filename, snpTemp_filename, 32);
	free(snp_kmers);

	std::cerr << "\n** K-mer create, now I write them on the binary file... \n" << std::endl;
}

/**
 * @info
 * @f1
 * @f2
 */
bool ExtractFrequency(std::string info, float& f1, float& f2)
{
	size_t startPos = info.find("CAF");
	char current_character;
	size_t i = 0;
	bool find_comma = false;

	if (startPos != std::string::npos)
	{
		//salto CAF=
		startPos += 4;

		size_t lastPos = startPos;

		while (!find_comma)
		{
			current_character = info[lastPos];
			if (current_character == ',')
				find_comma = true;
			else
			{
				lastPos++;
				i++;
			}
		}

	}
	else
		return false;

	//Get the first frequency
	std::string first_frequency = info.substr(startPos, i);

	f1 = std::stof(first_frequency);
	f2 = 1 - f1;
	return true;
}

/**
 *@snp_filename
 */
void convert_ucst_to_vcf(std::string snp_filename)
{
	//We open the stream to read the list of SNPs inside .txt file
	std::ifstream snp_file(snp_filename);

	//Create a new file whit extension .vcf in the same path were is present file whit snp in .txt format
	std::vector < std::string > columns = split_vargeno(snp_filename, '.');
	std::string snp_extension = columns[columns.size() - 1];
	std::string path_vcf = columns[0] + ".vcf";

	std::ofstream vcfFormat(path_vcf);

	//Check if the flow is active
	if (!snp_file.good())
	{
		std::cerr << "Error opening: " << snp_filename << std::endl;
		return -1;
	}

	if (!vcfFormat.good())
	{
		std::cerr << "Error create: " << path_vcf << std::endl;
		exit (EXIT_FAILURE);
	}

#define CHROM_FIELD   1//Name chromosome
#define INDEX_FIELD   2//Position within the CHROMOSOME where the SNP is located
#define STRAND_FIELD  6
#define REF1_FIELD    7//Base di riferimento di refNCBI
#define REF2_FIELD    8//Base di riferimento di refUCSC(dovrebbero essere uguali)
#define ALT_FIELD     9//Campo del tipo A/G in cui una base è il riferimento mentre l'altra è la base alternativa
#define TYPE_FIELD    11
#define COUNT_FIELD   21
#define ALLELES_FIELD 22
#define FREQS_FIELD   24

	//Header file
	vcfFormat << "#CHROM \t" << "POS \t" << "ID \t" << "REF \t" << "ALT \t" << "QUAL \t" << "FILTER \t" << "INFO \t" << std::endl;

	//Record extracted from the .txt file
	std::string line_string;

	//Split line_string
	std::vector < std::string > line_split;

	std::string stringOut = "";

	//I analyze each line inside .txt file
	while (std::getline(snp_file, line_string))
	{
		/*
		 * Discard all the records related to Meta-information
		 * (for example fileformat,fileDate, reference ..)
		 */
		if (line_string.empty())
		{
			stringOut = "";
			continue;
		}
		if (line_string[0] == '#')
		{
			stringOut = "";
			continue;
		}

		line_split = split_vargeno(line_string, '\t');

		// Copy chromosome name into an independent buffer
		std::string chr_name = line_split[CHROM_FIELD];

		//Add
		stringOut += chr_name + '\t';

		//Position within the CHROMOSOME where the SNP is located
		int unsigned index = std::stoi(line_split[INDEX_FIELD]);  // 0-based

		index++; //For .vcf is correct 1-based

		//Add position
		stringOut += std::to_string(index) + '\t';

		//Add ID
		stringOut += ". \t";

		//Read REF base
		const char ref_base = toupper(line_split[REF1_FIELD][0]);

		const bool neg = (line_split[STRAND_FIELD][0] == '-');
		if (!neg)
		{
			if (line_split[STRAND_FIELD][0] != '+')
			{
				std::cerr << "Error snp-file" << std::endl;
				exit (EXIT_FAILURE);
			}
		}

		const char a1 = neg ? rev(toupper(line_split[ALLELES_FIELD][0])) : toupper(line_split[ALLELES_FIELD][0]);
		const char a2 = neg ? rev(toupper(line_split[ALLELES_FIELD][2])) : toupper(line_split[ALLELES_FIELD][2]);
		if (!(a1 == 'A' || a1 == 'C' || a1 == 'G' || a1 == 'T'))
		{
			stringOut = "";
			continue;
		}
		if (!(a2 == 'A' || a2 == 'C' || a2 == 'G' || a2 == 'T'))
		{
			stringOut = "";
			continue;
		}
		if (a1 != ref_base && a2 != ref_base)
		{
			stringOut = "";
			continue;
		}

		stringOut += a1;
		stringOut += '\t';

		stringOut += a2;
		stringOut += '\t';

		std::string field_frequency = line_split[FREQS_FIELD];
		std::vector < std::string > v = split_vargeno(field_frequency, ',');

		float freq1 = std::stof(v[0]);
		float freq2 = std::stof(v[1]);
		/*
		 * freq1 non sempre corrisponde alla base di riferimento.
		 * Questo perchè possiamo avere(come nella prima righa del sito)
		 * refNCBI --- alleles      alleleFreqs
		 *    G		    A,G	     0.144169,0.855831
		 *
		 * freq1 non corrisponde alla frequenza di G e quindi
		 * devo fare lo swap
		 */
		if (a2 == ref_base)
		{
			const float tmp = freq1;
			freq1 = freq2;
			freq2 = tmp;
		}

		//Add QUAL
		stringOut += ". \t";

		//Add FILTER
		stringOut += ". \t";

		//Add INFO
		stringOut = stringOut + "CAF=" + std::to_string(freq1) + "," + std::to_string(freq2) + '\t';
		vcfFormat << stringOut << std::endl;
		stringOut = "";
	}

#undef CHROM_FIELD
#undef INDEX_FIELD
#undef STRAND_FIELD
#undef REF1_FIELD
#undef REF2_FIELD
#undef ALT_FIELD
#undef TYPE_FIELD
#undef COUNT_FIELD
#undef ALLELES_FIELD
#undef FREQS_FIELD

	snp_filename = path_vcf;
	vcfFormat.flush();
	vcfFormat.close();
}
