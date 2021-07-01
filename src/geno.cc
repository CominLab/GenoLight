#include <iostream>
#include <assert.h>
#include <vector>
#include <limits>
#include <unordered_map>
#include <string>
#include <cstring>
#include "geno.h"
#include "utility.h"
#include "mapping_reads.h"

int verbose_geno = 0;

param_creation_geno geno_param_init(int argc, char *argv[]);
std::vector<chrlens*> load_charlens_file(std::string chrlens_filename);
static void create_genotype_output_file(param_creation_geno* geno_param, unordered_map<string, pair<char, double>>* snp_2_genotype);

/**
 *@geno_param
 *@pileup_table
 *@pileup_size
 */
void genotype(param_creation_geno* geno_param)
{
	//uint32_t max_pos = rsd->max_pos;
	//uint32_t ambig_hits = 0;

	//Load chrlens file
	std::vector<chrlens*> chrlens = load_charlens_file(geno_param->chrlens_filename);

	/*
	 * The key is a string with the following format: name_chromosome + "$" + index_snp
	 * The value is a pair of <char,double>. The first can take as a value {0,1,2}
	 * 0-> REF/REF
	 * 1->REF/ALT
	 * 2->ALT/ALT
	 *
	 * while the second is the confidence of the call.
	 */
	unordered_map<string, pair<char, double>>* snp_2_genotype = new unordered_map<string, pair<char, double>>();

	size_t ref_call_count = 0;
	size_t alt_call_count = 0;
	size_t het_call_count = 0;

	struct pileup_entry **table = (geno_param->ptable).table;
	const size_t pileup_size = (geno_param->ptable).size;

	//for (size_t i = 0; i < geno_param->pileup_size; i++)
	for (size_t i = 0; i < pileup_size; i++)
	{
		for (struct pileup_entry *p = table[i]; p != NULL; p = p->next)
		{

			//struct packed_pileup_entry *p = &geno_param->pileup_table[i];

			//If the alternative base and the reference base are the same, it means that there is no SNP in that location
			if (p->ref == p->alt)
				continue;  // no SNP here

			//Due to the 'continue' we arrive here only in correspondence of the snp

			/*
			 * Warning: Initially variable "index" represents the position of the snp within the reference genome,
			 * 			then it contains the position of the snp within the current chromosome.
			 * 			For find current chromosome, I have to use the chrlens vector: in this
			 * 			way I can subtract the lengths of the previous chromosomes from the variable 'i'.
			 *
			 */

			size_t index = p->key;
			//size_t index = i;

			/*
			 * Variable j contains the index within the chrlens vector relative to the chromosome where
			 * the current snp is present
			 */
			size_t j;
			for (j = 0; j < chrlens.size() && index > chrlens[j]->len; j++)
				index -= chrlens[j]->len;

			/*
			 * Struct call have this member:
			 *  -genotype
			 *  -confidence
			 */
			const struct call call = choose_best_genotype(p->ref_cnt, p->alt_cnt, p->ref_freq, p->alt_freq);

			//unordered_map 'snp_2_genotype' key
			string snp_index = string(chrlens[j]->name) + "$" + to_string(index);

			//unordered_map 'snp_2_genotype' value
			pair<char, double> genotype_pair;

			switch (call.genotype)
			{
				case GTYPE_NONE:
					break;

				case GTYPE_REF:
					++ref_call_count;
					//std::cerr << "(105)Add snp_index: " << snp_index << std::endl;
					genotype_pair = make_pair('0', call.confidence);
					(*snp_2_genotype)[snp_index] = genotype_pair;
					break;

				case GTYPE_ALT:
					++alt_call_count;
					//std::cerr << "(107)" << std::endl;
					genotype_pair = make_pair('2', call.confidence);
					(*snp_2_genotype)[snp_index] = genotype_pair;
					break;

				case GTYPE_HET:
					++het_call_count;
					//std::cerr << "(114)" << std::endl;
					genotype_pair = make_pair('1', call.confidence);
					(*snp_2_genotype)[snp_index] = genotype_pair;
					break;
			}
		}
	}

	create_genotype_output_file(geno_param, snp_2_genotype);

	for (int i = 0; i < chrlens.size(); i++)
		delete chrlens[i];

	delete snp_2_genotype;
}

/**
 * Method in which the output .vcf file is created with the results of genotyping.
 *
 * @geno_param: Reference to a struct that contains all information concerning the genotyping phase
 * @snp_2_genotype: Reference to unordered_map in which the key is a string whit the format
 * 					[name_chromosome + "$" + index snp]. The return value is a pair of <char,double>.
 * 					The first can take as a value{0,1,2} [0->REF/REF, 1 ->REF/ALT, 2->ALT/ALT] while the
 * 					second is the confidence of the call.
 */
static void create_genotype_output_file(param_creation_geno* geno_param, unordered_map<string, pair<char, double>>* snp_2_genotype)
{
	//I open a channel with input file contains snp
	std::ifstream input(geno_param->path_snp);
	if (!input.good())
	{
		std::cerr << "Error opening: " << geno_param->path_snp << " . You have failed." << std::endl;
		return;
	}

	//I open a channel with output file
	std::ofstream output;
	output.open(geno_param->path_output);

	if (!output.good())
	{
		std::cerr << "Error opening: " << geno_param->path_output << " . You have failed." << std::endl;
		return;
	}

	//As we will see, iteratively contains all the lines of the input file or of the file containing the snp passed in input
	string line;

	/*
	 * Next two variables has_gt and  has_gq define if, inside .vcf file, in metadata section,
	 * there are the lines ##FORMAT=<ID=GT..> and <##FORMAT=<ID=GQ ..>.
	 *
	 * If such lines are present, I make them in the output file along with all the other metadata.
	 * Otherwise, as we will see, I will add them.
	 *
	 * Of course, the presence of these lines implies that there is also the FORMAT column (like # CHROM, POS, ID ...)
	 */
	bool has_gt = false;
	bool has_gq = false;

	/*
	 * Se all'interno nel file .vcf sono presenti tali righe, queste due variabili
	 * assumono il valore -1, altrimenti,come vedremo, gt_index = 0 e gq_index = -1
	 */
	int gt_index = -1;
	int gq_index = -1;
	bool head_has_gt_col = true;

	/**
	 * Iteratively, I consider all the lines in the input file(whit format .vcf or .txt if uscd)
	 * Basically the header lines are simply copied into the output file. For each SNP inside .vcf o .txt,
	 * check if, with input reads, it has been detected or not. This information is inside snp_2_genotype.
	 */
	while (std::getline(input, line))
	{
		if (line.empty())
			continue;

		if (line[0] == '#' and line[1] == '#')
		{
			//I copy the metadata from the .vcf file to genotype.out
			output << line << endl;

			/**
			 * As mentioned above, in the case where the lines ##FORMAT = <ID = GT, ..> and
			 * ##FORMAT = <ID = GQ, ..> are encountered set the relative variable at true
			 */
			if (line.find("ID=GT,") != std::string::npos)
				has_gt = true;

			else if (line.find("ID=GQ,") != std::string::npos)
				has_gq = true;

			continue;
		}
		else if (line[0] == '#')
		{
			/*
			 * I've finished the metadata section inside the file .vcf and line contain the header row
			 * (#CHROM,POS,ID,REF..)
			 */

			//As mentioned above, if these lines are not present I will add them
			if (!has_gt)
			{
				//Add this row inside output
				output << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;

				// randomly set >= 0
				gt_index = 0;
			}
			if (!has_gq)
			{
				//Add this row inside output
				output << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << endl;

				// randomly set >= 0
				gq_index = 1;
			}

			//output << line << endl;

			//Each element of the head_colimns vector contains an element of the snp header
			auto head_columns = split_vargeno(line, '\t');

			if (head_columns.size() < 10)
			{
				head_has_gt_col = false;

				/*
				 * The constant value GT:GQ corresponds to the FORMAT column.
				 *
				 * At the column DONOR corresponds the value,for example, 1/1 848,
				 * where 1/1 it means that for that SNP, both alleles have in the
				 * corresponding locus ALT base(0/0 both REF base while 0/1 one REF
				 * one ALT).
				 */
				line += "\tFORMAT\tDONOR";
			}

			//Copy this line into output
			output << line << endl;
			continue;
		}
		//Finally, I begin to consider all the lines containing the actual data of the snp

		//Each element of the vector 'columns' contains a value of snp
		vector < string > columns = split_vargeno(line, '\t');
		string chr_name = columns[0];

		if (chr_name[0] != 'c')
			chr_name = "chr" + chr_name;

		//Get position snp
		string pos = columns[1];

		//I remember this is the key to access the unordered map snp_2_genotype
		string snp_index = chr_name + "$" + pos;

		//std::cerr << "snp_index: " << snp_index << std::endl;

		if ((*snp_2_genotype).find(snp_index) == (*snp_2_genotype).end())
		{
			/**
			 * This means that with the reads provided as input, the snp contained in
			 * the chromosome 'chr_name' at the position 'pos' was not detected
			 */
			continue;
		}

		//If found, get value
		pair<char, double> genotype_pair = (*snp_2_genotype)[snp_index];

		//Default is REF/REF
		string genotype_string = "0/0";

		if (genotype_pair.first == '1')
			genotype_string = "0/1";

		else if (genotype_pair.first == '2')
			genotype_string = "1/1";

		int genotype_quality = -1 * 10 * log(genotype_pair.second);

		string format_str = "";

		/*
		 * If the FORMAT column (8) is present, I access its contents and store it in variable string 'format_str',
		 * then access the contents of the column (9) and store its contents in the variable 'info_str'.
		 */
		if (head_has_gt_col)
			format_str = columns[8];

		string info_str = "";

		if (head_has_gt_col)
			info_str = columns[9];

		//Split string 'format_str' and 'info_str'
		vector < string > format_columns;
		if (head_has_gt_col)
			format_columns = split_vargeno(format_str, ':');

		vector < string > info_columns;

		if (head_has_gt_col)
			info_columns = split_vargeno(info_str, ':');

		if (gt_index == -1 && has_gt)
		{
			//Dovrebbe essere ridondante, solo se has_gt == true gt_index==-1..
			for (std::vector<int>::size_type i = 0; i < format_columns.size(); i++)
			{
				if (format_columns[i] == "GT")
				{
					//indice all'interno di format_columns di GT, il valore corrispondente a GT in
					//info_columns è sempre all'indice gt_index
					gt_index = i;
					break;
				}
			}
			assert(gt_index >= 0);
		}
		if (gt_index == -1 && has_gq)
		{
			for (std::vector<int>::size_type i = 0; i < format_columns.size(); i++)
			{
				if (format_columns[i] == "GQ")
				{
					gq_index = i;								//indice all'interno di format_column di GQ, il valore corrispondente a GQ in
																//info_columns è sempre all'indice gt_index
					break;
				}
			}
			assert(gq_index >= 0);
		}
		if (has_gt)
		{
			info_columns[gt_index] = genotype_string;								 //Sostituisco il valore vecchio con quello nuovo
		}
		else
		{
			//Altrimenti lo aggiungo
			format_columns.push_back("GT");								 //Aggiungo in quanto non c'è
			info_columns.push_back(genotype_string);
		}

		if (has_gq)
		{
			//I replace the old value with the new one
			info_columns[gq_index] = to_string(genotype_quality);
		}
		else
		{
			format_columns.push_back("GQ");				//Aggiungo in quanto non c'è
			info_columns.push_back(to_string(genotype_quality));
		}

		// Create new format_str and info_str, connect them into one

		/*
		 * This new format is,for example: GT:GQ	0/0:846
		 */
		string new_format = format_columns[0];
		for (std::vector<int>::size_type i = 1; i < format_columns.size(); i++)
		{
			new_format += ":" + format_columns[i];
		}
		string new_info = info_columns[0];
		for (std::vector<int>::size_type i = 1; i < info_columns.size(); i++)
		{
			new_info += ":" + info_columns[i];
		}

		string new_line = columns[0];
		if (head_has_gt_col)
		{
			//Se ci sono modifico
			columns[8] = new_format;
			columns[9] = new_info;
		}
		else
		{
			//Altrimenti le aggiungo
			columns.push_back(new_format);
			columns.push_back(new_info);
		}
		for (std::vector<int>::size_type i = 1; i < columns.size(); i++)
		{
			new_line += '\t' + columns[i];
		}

		//I add to the output file a new record related to the snp detected by the reads provided in input
		output << new_line << endl;
	}

	//All lines of the input snp file have been analyzed

	//I close the channels
	input.close();
	output.close();
}
/**
 * I extract the input parameters from the input command supplied by the user
 *
 */
param_creation_geno geno_param_init(int argc, char *argv[])
{
	param_creation_geno geno_param;
	geno_param.n_threads = 1;
	geno_param.chunk_size = 1000;
	geno_param.quality_threshold = CONST_QUALITY_SCORE + 8;
	geno_param.inner_distance = 0;
	geno_param.read_len = 101;
	geno_param.path_output = "genotype.out";

	std::string nameVcfFile;

	//read the parameters
	for (int i = 1; i < argc; ++i)
	{
		if (strcmp(argv[i], "-c") == 0)
		{
			std::cerr << "-c" << std::endl;
			//k-mer length
			++i;
			geno_param.chunk_size = strtoull(argv[i], nullptr, 10);
		}
		else if (strcmp(argv[i], "-q") == 0)
		{
			++i;
			size_t qual = strtoull(argv[i], nullptr, 10);
			if ((qual >= 0) && (qual <= 93))
				geno_param.quality_threshold = CONST_QUALITY_SCORE + strtoull(argv[i], nullptr, 10);
			else
			{
				std::cerr << "Quality value must be between [0,93]" << std::endl;
				exit (EXIT_FAILURE);
			}
		}
		else if (strcmp(argv[i], "-t") == 0)
		{
			++i;
			geno_param.n_threads = strtoull(argv[i], nullptr, 10);

			if (geno_param.n_threads == 0)
				geno_param.n_threads = 1;
		}
		else if (strcmp(argv[i], "-r") == 0)
		{
			++i;
			geno_param.path_fastq = argv[i];
		}
		else if (strcmp(argv[i], "-i") == 0)
		{
			++i;
			geno_param.inner_distance = strtoull(argv[i], nullptr, 10);
		}
		else if (strcmp(argv[i], "-rl") == 0)
		{
			++i;
			geno_param.path_fastq_left = argv[i];
		}
		else if (strcmp(argv[i], "-rr") == 0)
		{
			++i;
			geno_param.path_fastq_right = argv[i];
		}
		else if (strcmp(argv[i], "-s") == 0)
		{
			++i;
			geno_param.path_snp = argv[i];
			string s = get_file_from_path(geno_param.path_snp);
			vector < string > columns = split_vargeno(s, '.');
			string snp_extension = columns[columns.size() - 1];
			nameVcfFile = columns[columns.size() - 2];

			if (snp_extension == "txt")
			{
				std::cerr << "\n Please pass the converted snp file in the .vcf format" << std::endl;
			}
		}
		else if (strcmp(argv[i], "-o") == 0)
		{
			++i;
			geno_param.path_output = std::string(argv[i]) + "genotype.out";
		}
		else if (strcmp(argv[i], "-p") == 0)
		{
			++i;
			geno_param.prefix = argv[i];
		}
		else
			print_help_geno();
	}

	if ((geno_param.path_fastq_right == "") && (geno_param.path_fastq_left != ""))
	{
		std::cerr << "Error: Please enter the path of the right reads! "<< std::endl;
		exit(EXIT_FAILURE);
	}
	else if ((geno_param.path_fastq_right != "") && (geno_param.path_fastq_left == ""))
	{
		std::cerr << "Please enter the path of the left reads!" << std::endl;
		exit(EXIT_FAILURE);
	}

	else if (((geno_param.path_fastq_right != "") || (geno_param.path_fastq_left != "")) && (geno_param.inner_distance == 0))
	{
		std::cerr << "Please enter the inner distance!" << std::endl;
		exit(EXIT_FAILURE);
	}

	geno_param.chrlens_filename = "output/" + geno_param.prefix + ".chrlens";
	geno_param.snp_dict_filename = "output/" + geno_param.prefix + "SmartDict" + nameVcfFile;
	
	//std::cerr << geno_param.snp_dict_filename << std::endl;

	std::string path_index_reassembly_snp = "output/" + geno_param.prefix + "Reassembly" + nameVcfFile;
	std::string path_index_ref_genome = "output/" + geno_param.prefix;

	//std::cerr << path_index_ref_genome << std::endl;
	//std::cerr << path_index_reassembly_snp << std::endl;

	geno_param.indexReferenceGenome = lib_aln_idx_load(path_index_ref_genome.c_str());
	geno_param.indexSNPsReassembly = lib_aln_idx_load(path_index_reassembly_snp.c_str());

	return geno_param;
}

std::vector<chrlens*> load_charlens_file(std::string chrlens_filename)
{
	std::string line;
	ifstream chrlens_file(chrlens_filename);
	if (!chrlens_file.good())
	{
		std::cerr << "Error open: " << chrlens_filename << std::endl;
		print_help_geno();
		exit (EXIT_FAILURE);
	}

	//La prima linea contiene il numero di basi del genoma ma a me ciò adesso non interessa
	std::getline(chrlens_file, line);

	std::getline(chrlens_file, line);
	size_t num_chromosomes = std::stoi(line);

	std::vector<chrlens*> chrlens_v(num_chromosomes);

	for (std::vector<int>::size_type i = 0; i < num_chromosomes; i++)
	{
		std::getline(chrlens_file, line);
		int pos = line.find(" ");
		std::string name = line.substr(0, pos);
		size_t len = std::stoi(line.substr(pos + 1, line.size() - 1));
		chrlens_v[i] = new chrlens(name, len);
	}
	chrlens_file.close();
	return chrlens_v;
}

static inline struct call choose_best_genotype(const int ref_cnt, const int alt_cnt, const uint8_t ref_freq_enc, const uint8_t alt_freq_enc)
{
	/*
	 * G0: 'HOMOZYGOUS'
	 * G1: 'HETEROZYGOUS'
	 * G3: 'HOMOZYGOUS ALTERNATE'
	 */
	static struct
	{
			double g0;  // P(counts|G0)/(ref_cnt + alt_cnt choose ref_cnt)
			double g1;  // P(counts|G1)/(ref_cnt + alt_cnt choose ref_cnt)
			double g2;  // P(counts|G2)/(ref_cnt + alt_cnt choose ref_cnt)
	} cache[MAX_COV + 1][MAX_COV + 1];

	/*
	 * This is used to penalize SNPs that have abnormally high coverages
	 * as well as those have lower coverages than expected
	 */
	static double poisson[2 * MAX_COV + 1];

	static bool init = false;

	if (!init)
	{
		for (int ref_cnt = 0; ref_cnt <= MAX_COV; ref_cnt++)
		{
			for (int alt_cnt = 0; alt_cnt <= MAX_COV; alt_cnt++)
			{
				cache[ref_cnt][alt_cnt].g0 = pow(1.0 - ERR_RATE, ref_cnt) * pow(ERR_RATE, alt_cnt);
				cache[ref_cnt][alt_cnt].g1 = pow(0.5, ref_cnt + alt_cnt);
				cache[ref_cnt][alt_cnt].g2 = pow(ERR_RATE, ref_cnt) * pow(1.0 - ERR_RATE, alt_cnt);
			}
		}

		//AVG_COV is average coverage of reads( costante pari a 7.1)
		const double M = exp(-AVG_COV);

		/*
		 * Viene calcolato poisson per ogni possibile alfa+beta
		 *  exp(lgamma(i + 1.0)) è equivalente a (i)!
		 */
		for (int i = 0; i <= (2 * MAX_COV); i++)
		{
			poisson[i] = (M * pow(AVG_COV, i)) / exp(lgamma(i + 1.0));
		}

		init = true;
	}

	if ((ref_cnt == 0 && alt_cnt == 0) || (ref_cnt == MAX_COV && alt_cnt == MAX_COV))
	{
		return CALL(GTYPE_NONE, 0.0);
	}
	const double g0 = cache[ref_cnt][alt_cnt].g0;
	const double g1 = cache[ref_cnt][alt_cnt].g1;
	const double g2 = cache[ref_cnt][alt_cnt].g2;

	/*
	 	 * Remember: p frequency of reference allele
	 	 * 			 q frequency of alternative allele
	*/

	const double p = ref_freq_enc / 255.0;
	const double q = alt_freq_enc / 255.0;

	const double p2 = p * p;
	const double q2 = q * q;

	const double p_g0 = p2 * g0;
	const double p_g1 = (1.0 - p2 - q2) * g1;
	const double p_g2 = q2 * g2;
	const double total = p_g0 + p_g1 + p_g2;

	const int n = ref_cnt + alt_cnt;

	/*
	 * CALL ritorna una struct contenente il tipo di genotipizzazione
	 * (GTYPE_REF,GTYPE_HET,GTYPE_ALT) e la confidenza della chiamata.
	 */
	if (p_g0 > p_g1 && p_g0 > p_g2)
	{
		return CALL(GTYPE_REF, ((double)(p_g0/total))*poisson[n]);
	}
	else if (p_g1 > p_g0 && p_g1 > p_g2)
	{
		return CALL(GTYPE_HET, ((double)(p_g1/total))*poisson[n]);
	}
	else
	{
		return CALL(GTYPE_ALT, ((double)(p_g2/total))*poisson[n]);
	}
}
