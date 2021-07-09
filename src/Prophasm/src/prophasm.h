#include<string>

void print_help_prophasm();

uint64_t goProphasmReassembly(uint8_t argc, std::string path_file_kmer, std::string path_output);

void test_file(FILE *fo, std::string fn);

template<typename _nkmer_T> 
int32_t encode_forward(const char *kmers, const int32_t k, _nkmer_T &nkmer);

template<typename _nkmer_T>
int32_t encode_reverse(const char *kmers, const int32_t k, _nkmer_T &nkmer);

template<typename _nkmer_T>
int32_t encode_canonical(const char *kmers, const int32_t k, _nkmer_T &nkmer);

template<typename _nkmer_T>
int32_t decode_kmer(_nkmer_T nkmer, int32_t k, std::string &kmer);

void reverse_complement_in_place(std::string &kmer);

template<typename _set_T>
void debug_print_kmer_set(_set_T &set, int k, bool verbose);

template<typename _set_T>
void debug_print_kmer_set(_set_T &set, int k, bool verbose);

struct contig_t;
    
template<typename _set_T>
int kmers_from_fasta(const std::string &fasta_fn, _set_T &set, int32_t k, FILE* fstats,bool verbose);

template<typename _set_T>
int32_t find_intersection(const std::vector<_set_T> &sets, _set_T &intersection);

template<typename _set_T, typename _subset_T>
int32_t remove_subset(std::vector<_set_T> &sets, const _subset_T &subset);

template<typename _set_T>
uint64_t assemble(const std::string &fasta_fn, _set_T &set, int32_t k, FILE* fstats, bool verbose);
