#include <iostream>
#include <string>
#include "fq.hpp"
#include <limits>
#include <vector>
#include <fstream>

//costruttore
//https://stackoverflow.com/questions/50417435/why-cannot-i-initialize-a-ifstream-ref-from-constructor
fq::fq(const std::string &_path_fastq) :
fastq(_path_fastq)
{
	if (_path_fastq != " ")
		path_fastq = _path_fastq;
	else
	{
		std::cerr << "Please enter the path where the reads are stored." << std::endl;
		exit (EXIT_FAILURE);
	}
	//fastq(path_fastq); //std::ios::binary
}

//costruttore
fq::fq(std::string f_filename_fastq, std::string s_filename_fastq)
{
	if ((f_filename_fastq != " ") && (s_filename_fastq != " "))
	{
		path_fastq = f_filename_fastq;
		path_fastq_mate_pairs = s_filename_fastq;
	}
	else
	{
		std::cerr << "Please enter the path where the reads mate pairs stored." << std::endl;
		exit (EXIT_FAILURE);
	}
	std::ifstream fastq(path_fastq); //std::ios::binary
}

FastqRecord fq::read_fastq_record(std::istream& fastq)
{
	std::string header;
	std::string seq;
	std::string quals;

	// Unlike FASTA, FASTQ always use one line per field.
	if (!std::getline(fastq, header) || header == "")
	{
		end_file = true; //private variable
		FastqRecord f = FastqRecord();
		return f;
	}
	std::getline(fastq, seq);
	fastq.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	std::getline(fastq, quals);
	if (seq.size() != quals.size())
	{
		std::cerr << "Length read different length quality score" << std::endl;
		exit (EXIT_FAILURE);
	}

	FastqRecord f = FastqRecord(std::move(header), std::move(seq), std::move(quals));
	return f;
}

void fq::fill_chunk(std::vector <FastqRecord>& chunks)
{
	int i = 0;
	for (; ((i < chunks.size()) && (end_file == false)); i++)
	{
		chunks[i] = read_fastq_record(fastq);
	}
	if (end_file == false)
		chunks.resize(i);
	else
	{
		chunks.resize(i - 1);//If end_file==true add a nullptr in records
		//std::cerr << "Ultimo chunk size: "<< records.size() << std::endl;
		fastq.close();
	}
}

bool fq::AllReadExtracted()
{
	return end_file;
}
