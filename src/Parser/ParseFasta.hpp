#ifndef PARSE_FASTA_HPP
#define PARSE_FASTA_HPP

#include <map>
#include <string>
#include <vector>
#include <fstream>

struct FastaRecord
{
		std::string name;
		std::string sequence;

		FastaRecord();

		FastaRecord(std::string&& name, std::string&& sequence) : name {name}, sequence {sequence}
		{
		}
};

class ParseFasta
{
	public:
		//constructor
		ParseFasta(std::string filename_fasta);

		//METODI PUBLICI
		FastaRecord* find_seq_by_name(const std::string current_chr_name, unsigned int& start_index_chr);

		FastaRecord* GetSequence(size_t i);

		bool IsChromosomePresent(std::string chromosome_name);

		uint64_t GetNumOfGenomeBases();

		size_t Size();

		//~ParseFasta();

	private:
		static constexpr char fasta_delim = {'>'};

		//STRUTTURA UTILIZZATA PER MEMORIZZARE COPPIA NOME CROMOSOMA-SEQUENZA
		std::vector<FastaRecord*> chromosome {};

		std::string filename_fasta;

		uint64_t number_of_genome_bases = 0;

		//PRIVATE METHOD

		std::size_t count_records(std::istream& fasta);

		void ExtractNameFromDescription(std::string& name);

		FastaRecord* read_fasta_record(std::istream& fasta);

		size_t count_records(std::istream& is, const char record_delim);

		void read_fasta(const std::string& fasta_path);
};

#endif
