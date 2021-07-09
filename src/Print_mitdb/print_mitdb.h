#ifndef PRINTMITDB_H_
#define PRINTMITDB_H_

#define POW_2_32 (1UL << 32)

#include <string.h>
#include <fstream>
#include <iostream>

using namespace std;

typedef u_int64_t pac_t;
typedef u_int32_t pos_t;
typedef u_int8_t flag_t;

int read_vargeno_ref_dict (string path_to_vargeno_ref, string path_to_out, string path_smartDictionary);
int read_vargeno_snp_dict(string path_to_snp_dict, string path_to_out, string path_smartDictionary);
string decode_kmer_from_binary(pac_t e);

#endif
