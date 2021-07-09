#ifndef FASTQ_HPP
#define FASTQ_HPP

#include <map>
#include <string>
#include <vector>
#include <fstream>

struct FastqRecord
{
		std::string header;
		std::string seq;
		std::string qua;

		FastqRecord()
		{
			header = " ";
			seq = " ";
			qua = " ";
		}

		FastqRecord(std::string _header, std::string _seq, std::string _qua)
		{
			header = _header;
			seq = _seq;
			qua = _qua;
		}
};

class fq
{
	public:
		//constructor
		fq(const std::string &filename_fastq);
		fq(std::string f_filename_fastq, std::string s_filename_fastq);

		//METODI PUBLICI
		void fill_chunk(std::vector<FastqRecord>& chunk);
		bool AllReadExtracted(void);
		//~ParseFasta();

	private:
		static constexpr char fastq_delim = {'@'};
		uint64_t index_last_read;
		std::string path_fastq;
		std::string path_fastq_mate_pairs;
		bool end_file = false;
		std::ifstream fastq;
		//PRIVATE METHOD
		FastqRecord read_fastq_record(std::istream& fastq);
};

#endif
