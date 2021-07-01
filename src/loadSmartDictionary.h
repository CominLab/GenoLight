#ifndef LOAD_SMART_DICTIONARY
#define LOAD_SMART_DICTIONARY

#include <iostream>
#include <fstream>
#include <string>

#include "utility.h"
#include "dictgen.h"
#include "genolight.h"
#include "pileup.h"

using namespace std;

//Definiscono la dimension,di un byte, di un record memorizzato all'interno dell dizionario nel file .bin
#define BYTE_REC_REF 5
#define BYTE_STORE_REC_REF 9
#define BYTE_REC_SNP 13


#define BYTE_HEADER_SMART_DICT 33
#define BYTE_STORE_POS 4

#define POSITION_NUM_COMPLEATE_KMER 25

#define COMPLETE_SNP_DICT 12 //(1+1+1+1 + 4 + 4)

SNPSmartDictionary* load_smartSNP_Dict(param_creation_geno* geno_param);
void complete_smart_dict(std::string path_smartDict, std::string path_reassembly_snp, std::string snpTemp_filename);
#endif

