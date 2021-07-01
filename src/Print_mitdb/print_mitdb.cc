#include <cstdlib>
#include <string>
#include <assert.h>

#include <iostream>
#include <fstream>

#include "print_mitdb.h"
#include "utility.h"

using namespace std;

int read_vargeno_ref_dict(string path_to_vargeno_ref, string path_to_output, string path_smartDictionary)
{
    pac_t dict_siz, aux_table_size;
    pac_t packed;

    pos_t pos_kmer_rif;

    //Lascio semplicemente il posto
    pos_t pos_kmer_reassembly = 0;

    flag_t amb_flag;
    pac_t i;
    string kmer;

    string decode_kmer;

    //Serve per leggere il contenuto del file binario
    ifstream fp (path_to_vargeno_ref, ios:: in | ios::binary);

    //Serve per scrivere il contenuto del file binario precedentemente letto
    ofstream fStoreKmer (path_to_output , ios:: out);

    //Questo serve per memorizzare le informazzioni come flag e pos della corrispondente k-mer
    ofstream fStoreInfoKmer (path_smartDictionary, ios::out | ios::binary);

    if(fp.is_open() == false)
    {
        cerr << "Unable to open the specified ref file" << path_to_vargeno_ref << endl;
        return 1;
    }   
    if(fStoreKmer.is_open() == false)
    {
           cerr << "Unable to create the specified ref file" <<  path_to_output << endl;
           return 2;
    }

    if(!fp.read((unsigned char*) &dict_siz, sizeof(pac_t))) //read the two uint64_t at the beginning (size and placeholder)
    {
        cerr << "Unable to read the total number of kmers for" << path_to_vargeno_ref << endl;
        return 3;
    }
    if(!fp.read((unsigned char*) &aux_table_size, sizeof(pac_t)))
    {
        cerr << "Unable to read the size of the auxiliary table for" << path_to_vargeno_ref << endl;
        return 4;
    }
    if (dict_siz > POW_2_32)
    {
		cerr << "Reference dictionary is too large (limit:" << POW_2_32 << "32-mers" << endl;
		return 5;
    }
    //CONTROLLO ANCHE SU QUELLA DI SUPPORTO
    
    for(i = 0; i < dict_siz; ++i)
    {
        packed = read_uint64_t(fp);
        decode_kmer = decode_kmer_from_binary(packed);
        fStoreKmer << ">k\n" << decode_kmer << endl; //Necessario mettere >k altrimenti non funziona,prophasm non riesce a leggere le k-mer

        pos_kmer_rif = read_uint32_t(fp);
        amb_flag = read_uint8_t(fp); //AMBIGUOUS FLAG

        fStoreInfoKmer.write((char *) &pos_kmer_reassembly, sizeof(pos_kmer_reassembly));
        fStoreInfoKmer.write((char *) &pos_kmer_rif, sizeof(pos_kmer_rif));
        fStoreInfoKmer.write((char *) &amb_flag, sizeof(amb_flag));
    }
    

    //Ricopio tabella ausiliaria
    for(i = 0; i < aux_table_size*10; ++i)
    {
    	pos_kmer_rif = read_uint32_t(fp);
    	fStoreInfoKmer.write((char *) &pos_kmer_reassembly, sizeof(pos_kmer_reassembly));
    }

    fp.close();
    fStoreKmer.close();
    fStoreInfoKmer.close();

    return 0;
}

int read_vargeno_snp_dict(string path_to_snp_dict, string path_to_output, string path_smartDictionary)
{
    pac_t dict_siz, aux_table_size;
    pac_t packed;

    pos_t pos_kmer_rif;
    pos_t pos_kmer_reassembly = 0;

    flag_t snp_info;
    flag_t amb_flag;

    flag_t ref_freq;
    flag_t alt_freq;
    
    pac_t i;
    string decode_kmer;

    
    ifstream fp (path_to_snp_dict, ios::in | ios::binary);
    ofstream fStoreKmer (path_to_output, ios:: out);

    ofstream fStoreInfoKmer (path_smartDictionary, ios::out | ios::binary);

    if(fp.is_open() == false)
    {
    	cerr << "Unable to open the specified SNP file" << path_to_snp_dict << endl;
        return 1;
    }   
    
    if(!fp.read((unsigned char*) &dict_siz, sizeof(pac_t))) //read the two uint64_t at the beginning (size and placeholder)
    {
        cerr << "Unable to read the total number of kmers for" << path_to_snp_dict << endl;
        return 2;
    }
    if(!fp.read((unsigned char*) &aux_table_size, sizeof(pac_t)))
    {
        cerr << "Unable to read the size of the auxiliary table for" << path_to_snp_dict << endl;
        return 3;
    }

    for(i = 0; i < dict_siz; ++i)
    {
    	//BUFFER DIVENTERA' UN VETTORE E SI COPIERA' QUELLO SU DISCO

        packed = read_uint64_t(fp);
        pos_kmer_rif = read_uint32_t(fp);
        snp_info = read_uint8_t(fp); //SNP info
        amb_flag = read_uint8_t(fp); //AMBIGUOUS FLAG
        ref_freq = read_uint8_t(fp); //ref freq
        alt_freq = read_uint8_t(fp); //alt freq

        decode_kmer = decode_kmer_from_binary(packed);

        fStoreKmer << ">k\n" << decode_kmer << endl; //Necessario mettere >k altrimenti non funziona,prophasm non riesce a leggere le k-mer"

        fStoreInfoKmer.write((char *)&pos_kmer_reassembly, sizeof(pos_kmer_reassembly));
        fStoreInfoKmer.write((char *)&pos_kmer_rif, sizeof(pos_kmer_rif));
        fStoreInfoKmer.write((char *)&snp_info, sizeof(snp_info));
        fStoreInfoKmer.write((char *)&amb_flag, sizeof(amb_flag));
        fStoreInfoKmer.write((char *)&ref_freq, sizeof(ref_freq));
        fStoreInfoKmer.write((char *)&alt_freq, sizeof(alt_freq));
    }
    
    for(i = 0; i < aux_table_size* 8; ++i)
    {
       	pos_kmer_rif = read_uint32_t(fp);
       	fStoreInfoKmer.write((char *) &pos_kmer_rif, sizeof(pos_kmer_rif));
    }

    fp.close();
    fStoreKmer.close();
    fStoreInfoKmer.close();

    return 0;
}

string decode_kmer_from_binary(pac_t e)
{
	string decode_kmer;

    int i;
	for (i=0; i<32; ++i)
	{
		switch (e & 0x0000000000000003)
        {
			case 0:
				decode_kmer += "A";
				break;
			case 1:
				decode_kmer += "C";
				break;
			case 2:
				decode_kmer += "G";
				break;
			case 3:
				decode_kmer += "T";
				break;
		}
		e >>= 2;
	}

	return decode_kmer;
}










