#include <cstdio>
#include <time.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <errno.h>
#include <sys/stat.h>

#include "dictgen.h"
#include "utility.h"

#include "lib_aln_inexact_matching.h"
#include "ParseFasta.hpp"
#include "fq.hpp"
#include "print_mitdb.h"
#include "prophasm.h"
#include "loadSmartDictionary.h"
#include "geno.h"
#include "mapping_reads.h"
#include "pileup.h"

#define PILEUP_TABLE_INIT_SIZE (1 << 25)

#define NUM_MANDATORY_PARAM_CREATE_INCOMPLEATE_DICTIONARY 2
#define NUM_MANDATORY_PARAM_REASSEMBLY 2
#define NUM_MANDATORY_PARAM_CREATE_FMDINDEX 1
#define NUM_MANDATORY_PARAM_COMPLETE_SMART_DICTIONARY 1
#define NUM_MANDATORY_PARAM_GENO 6
#define NUM_MANDATORY_PARAM_INDEX_VCF 10
#define NUM_MANDATORY_PARAM_INDEX_UCSC 9

void createIndex(std::string path)
{
	string file = get_file_from_path(path);
	vector < string > nameFileWhitExtension = split_vargeno(file, '.');
	string nameFile = "output/" + nameFileWhitExtension[0];
	lib_aln_index(path.c_str(), nameFile.c_str(), BWTALGO_AUTO);
}

int main(const int argc, const char *argv[])
{
	//Instant start of command execution
	clock_t begin;

	//Instant end of command execution
	clock_t end;

	//Variable defined as: end-begin
	double time_spent;

	/**
	 *
	 * The possible commands that the user can provide are:
	 *
	 * create_incomplete_smartSnpDictionary: Create an incomplete smart dictionary containing k-mer extracted from a reference genome in which
	 * 										  of snp which replaces the alternative allele at the base of reference
	 *
	 * reassembly: Command that starts the re-assembly of the k-mer dataset provided as input
	 *
	 * createFMDIndex: Command to perform the indexing of a string whose characters belong to the alphabet {A, C, G, T}
	 *
	 * complete_smartSnpDictionary: Command that inserts, inside the position_to_reassembly field, the position of it inside the re-assembly.
	 *
	 * geno: Command to perform genotyping
	 *
	 *
	 */

	if (argc < 2)
	{
		std::cerr << "\nCommand not inserted, please enter a command." << std::endl;

		//Show the set of valid commands to the user and interrupt the process
		print_help_valid_comand();
	}

	const std::string inputUserComand = argv[1];

	/* Command that produce the dictionary containing the k-mers extracted from the reference genome
	 * where replacing the reference base with the alternative base.
	 */
	if (inputUserComand.compare("create_incomplete_smartSnpDictionary") == 0)
	{
		//Path to reference genome
		string ref_filename;

		//Path of list of SNPs
		string snp_filename;

		std::string snp_dict_filename;

		const std::string path_output = "output/";

		//Path where store file in which to store the k-mer in the literal format
		string temp_filename;

		//Path in which to store file that contains k-mer in literal format
		string chrlens_filename;

		string prefix;
		if ((argc-1)<6)
		{
			print_help_incDic();
		}
		for (int i = 2; i < argc; ++i)
		{
			if (strcmp(argv[i], "-r") == 0)
			{
				++i;
				ref_filename = argv[i];
			}
			if (strcmp(argv[i], "-s") == 0)
			{
				++i;
				snp_filename = argv[i];
			}
			else if (strcmp(argv[i], "-p") == 0)
			{
				++i;
				prefix = argv[i];
			}
			else
			{
				std::cerr << "Parameter: " + std::string(argv[i]) + " not admissible." << std::endl;
				print_help_incDic();
			}
		}

		if (prefix.size() == 0)
			prefix = get_name_from_genome_path(ref_filename);

		create_default_out_folder(path_output);

		string s = get_file_from_path(snp_filename);
		vector < string > columns = split_vargeno(s, '.');
		s = columns[columns.size() - 2];

		temp_filename = path_output + prefix + "Temp.txt";
		chrlens_filename = path_output + prefix + ".chrlens";
		snp_dict_filename = path_output + prefix + "SmartDict" + s;

		cout << "\nStarting genome parsing..." << '\n';

		ParseFasta genome(ref_filename);

		cout << "Finished genome parsing!" << '\n';

		/*
		 * The head of the file contains the number of records present.
		 * Each record in the file .chrlen contains the name of the chromosome and its length.
		 * This file is useful when genotyping will be performed
		 *
		 */
		CreateChrlensFile(chrlens_filename, genome);

		cout << "\nStarting snp smart dictionary..." << '\n';
		/*
		 * This vector stores a Boolean value at the i-th index depending on whether or not a snp is at the i-th
		 * position of the reference genome.
		 * In terms of memory it is not the most efficient solution, however it allows in constant time to check if in a
		 * certain position there is a snp or not.
		 */
		uint64_t num_record_dictionary = 0;

		bool *snp_locations;
		size_t snp_locs_size = genome.GetNumOfGenomeBases();

		make_snp_dict_from_vcf(genome, snp_dict_filename, snp_filename, temp_filename, &snp_locations, snp_locs_size);

		cout << "Finished snp smart dictionary!" << '\n';
		//I create a file containing the value of k and read len with which the dictionaries were created
		//CreateInfoFile(path_info_file, len_reass_snp);
	}
	else if (inputUserComand.compare("reassembly") == 0)
	{
		//Path to store the file containing the reassembly
		string path_store_reassembly_snp;

		//
		string prefix;

		//
		string nameVcfFile;

		//read the parameters
		if ((argc-1)<4)
		{
			print_help_reassembly();
		}
		for (int i = 2; i < argc; ++i)
		{
			if (strcmp(argv[i], "-n") == 0)
			{
				++i;
				nameVcfFile = argv[i];
			}
			else if (strcmp(argv[i], "-p") == 0)
			{
				++i;
				prefix = argv[i];
			}
			else
			{
				std::cerr << "Parameter: " + std::string(argv[i]) + " not admissible." << std::endl;
				print_help_reassembly();
			}
		}

		std::cerr << "\nStarting re-assembly dictionary..." << std::endl;

		//Per default
		string snpTemp_filename = "output/" + prefix + "Temp.txt";

		//
		string path_smartDict = "output/" + prefix + "SmartDict" + nameVcfFile;

		//Devo aggiungere una estensione per simmetria
		std::string path_reass = "output/" + prefix + "Reassembly" + nameVcfFile + ".txt";

		//I use the Prophasm software to re-assemble k-mer inside snpTemp_filename
		uint64_t len_reass_snp = goProphasmReassembly(32, snpTemp_filename, path_reass);

		std::cerr << "Completed the re-assembly of dictionary. Length of re-assembly is: " << len_reass_snp << std::endl;

		//Used to write the re-assembly position inside smart dictionary and to update field flag inside dictionary
		ofstream smartDictionary_write(path_smartDict, ios::binary | ios::out | ios::in);

		std::cerr << "Length of re-assembly is: " << len_reass_snp << std::endl;

		//The first 16 bytes of the dictionary store the size of the dictionary and auxiliary table.
		//In addition, another byte is used to store the value of k.
		smartDictionary_write.seekp(17, ios_base::beg);
		smartDictionary_write.write((char *) (&len_reass_snp), sizeof(uint64_t));

		smartDictionary_write.flush();
		smartDictionary_write.close();
	}
	else if (inputUserComand.compare("createFMDIndex") == 0)
	{
		/*
		 std::cerr << "Indexing started.." << std::endl;
		 string s = get_file_from_path(argv[2]);
		 vector < string > columns = split_vargeno(s, '.');
		 s = "output/" + columns[columns.size() - 2];

		 lib_aln_index(argv[2], s.c_str(), BWTALGO_AUTO);
		 std::cerr << "Indexing end!" << std::endl;
		 */
		std::cerr << "Indexing..." << std::endl;
		createIndex(argv[2]);
		std::cerr << "Finished!" << std::endl;
	}
	else if (inputUserComand.compare("complete_smartSnpDictionary") == 0)
	{
		std::string prefix;

		//
		string nameVcfFile;

		//read the parameters
		if ((argc-1)<4)
		{
			print_help_compDic();
		}
		for (int i = 2; i < argc; ++i)
		{
			if (strcmp(argv[i], "-p") == 0)
			{
				++i;
				prefix = argv[i];
			}
			else if (strcmp(argv[i], "-n") == 0)
			{
				++i;
				nameVcfFile = argv[i];
			}
			else
			{
				std::cerr << "Parameter: " + std::string(argv[i]) + " not admissible." << std::endl;
				print_help_compDic();
			}
		}

		//Path where the dictionary is stored
		std::string snpDictionary_filename = "output/" + prefix + "SmartDict" + nameVcfFile;

		//Path containing the k-mer in the literal format
		std::string snpTemp_filename = "output/" + prefix + "Temp" + ".txt";

		//Path containing the re-assembly indexation
		std::string path_index_reassembly_snp = "output/" + prefix + "Reassembly" + nameVcfFile;

		std::cerr << "\nStarting complete SNP dictionary..." << std::endl;

		//Complete smart dictionary
		complete_smart_dict(snpDictionary_filename, path_index_reassembly_snp, snpTemp_filename);

		std::cerr << "\n Smart SNP dictionary completed successfully!\n" << std::endl;

		//std::remove(dict_param.snpTemp_filename.c_str());

		//I create a file containing the value of k and read len with which the dictionaries were created
	}
	else if (inputUserComand.compare("geno") == 0)
	{

		/*
		 * It is the last phase of the program. We have already created (and completed) the smart dictionary,
		 * created and indexed the re-assembly of the SNPs smart dictionary and the reference smart
		 * dictionary.
		 */

		//I extract the input parameters from the command supplied by the user regarding the genotyping phase
		param_creation_geno geno_param = geno_param_init(argc, argv);

		/*
		 * We assume that the maximum position encountered in the reference dictionary will not be smaller than the
		 * maximum position that will be encountered in the SNP dictionary. If not, we will reallocate our pileup table.
		 */
		geno_param.pileup_size = GetNumBasesReferenceGenome(geno_param.chrlens_filename);

		//std::cerr << "Size pileup table: " << (unsigned) geno_param.pileup_size << std::endl;

		fprintf(stderr, "\nINPUT PARAMETERS\n");
		fprintf(stderr, "Quality threshold: %zu \n", geno_param.quality_threshold - CONST_QUALITY_SCORE);
		fprintf(stderr, "Len reads: %zu \n", geno_param.read_len);

		//std::cerr << "\nInput parameters:" << std::endl;
		//std::cerr << "Quality threshold: " << (unsigned)geno_param.quality_threshold - CONST_QUALITY_SCORE << std::endl;
		//std::cerr << "Len reads: " << (unsigned)geno_param.read_len << "\n"<< std::endl;

		/*
		 if (geno_param.inner_distance == 0)
		 std::cerr << "Use single end reads: " << std::endl;
		 else
		 std::cerr << "Use mate pairs reads. Inner distance: " << (unsigned)geno_param.inner_distance << std::endl;
		 */
		//geno_param.pileup_table = (struct packed_pileup_entry*) calloc(geno_param.pileup_size, sizeof(*geno_param.pileup_table));
		ptable_init(&geno_param.ptable, PILEUP_TABLE_INIT_SIZE);

		std::cerr << "Loading smart snp dictionary...\n" << std::endl;
		SNPSmartDictionary* snpRD = load_smartSNP_Dict(&geno_param);
		std::cerr << "Finished loading!\n" << std::endl;

		//Using the reads, update the pilleup table
		begin = clock();

		std::cerr << "Starting reads mapping...\n" << std::endl;

		const uint8_t num_mismatches = 1;
		def_param* dp = lib_aln_set_default_parameters(geno_param.indexReferenceGenome, geno_param.indexSNPsReassembly, snpRD,
														geno_param.quality_threshold, num_mismatches);
		updatePilleupSingleEnd(&geno_param, dp);

		end = clock();
		time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
		printf("Time: %f sec\n", time_spent);

		free(dp);
		free(snpRD->snp_dict);

		//principalDictionaty_dealloc(snpRD->snp_dict);
		//free(snpRD->snp_dict);

		delete[] snpRD->snp_auxiliary_table;

		total_kmer_snp_struct_dealloc(snpRD->reverse_total_kmer_snp);
		free(snpRD->reverse_total_kmer_snp);

		delete snpRD;

		lib_aln_idx_destroy(geno_param.indexReferenceGenome);
		lib_aln_idx_destroy(geno_param.indexSNPsReassembly);

		//Starting the genotyping
		begin = clock();
		std::cerr << "Starting genotyping...\n" << std::endl;
		genotype(&geno_param);
		std::cerr << "Finished genotyping!\n" << std::endl;

		ptable_dealloc(&geno_param.ptable);
		//free(geno_param.ptable);

		end = clock();
		time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
		printf("Time: %f sec\n", time_spent);
		std::cout << "Memory peak (MB): " << to_string(getPeakVirtualMemoryUsed() / 1024) << std::endl;
	}
	else if (inputUserComand.compare("help") == 0)
	{
		print_help_valid_comand();
	}
}
