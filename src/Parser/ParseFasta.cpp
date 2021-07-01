#include <iostream>
#include <string>
#include "ParseFasta.hpp"
#include <limits>
#include <vector>
#include <fstream>

//Constructor
ParseFasta::ParseFasta(std::string _filename_fasta)
{
	if (filename_fasta != " ")
		filename_fasta = _filename_fasta;
	else
	{
		std::cerr << "Please enter the path where the reference genome is stored." << std::endl;
		exit (EXIT_FAILURE);
	}
	read_fasta(_filename_fasta);
}

//Destructor
/*
ParseFasta::~ParseFasta()
{
	for(int i = 0; i < chromosome.size(); i++)
	{
		delete[] chromosome[i];
	}
}
*/

// Counts the occurrences of record_delim that follow a newline '\n'
inline size_t ParseFasta::count_records(std::istream& is, const char record_delim)
{
	const auto current_position = is.tellg();
	size_t result {};
	while (is)
	{
		if (is.peek() == record_delim)
			++result;
		is.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	is.clear();
	is.seekg(current_position, std::ios::beg);
	return result;
}

/**
 *
 */
FastaRecord* ParseFasta::read_fasta_record(std::istream& fasta)
{
	std::string name;
	std::getline(fasta, name); // name is always a single line

	std::string line;
	std::getline(fasta, line);

	// The FASTA format is not as simple as FASTQ - the sequence
	// may be broken into multiple lines. We assume each line is the same size.
	if (!fasta.good() || fasta.peek() == fasta_delim)
	{
		return new FastaRecord(std::move(name), std::move(line));
	}
	else
	{
		const auto line_size = line.size();
		line.resize(line_size);

		std::string seq {};
		seq.reserve(line_size * 100);

		const auto line_begin = line.cbegin(), line_end = line.cend();
		seq.insert(seq.end(), line_begin, line_end);
		while (fasta.good())
		{
			fasta.getline(&line[0], line_size + 1);
			if (line.front() == fasta_delim)
			{
				if (fasta.good())
					fasta.seekg(-fasta.gcount(), std::ios_base::cur);
				break;
			}
			if (fasta.gcount() < line_size)
			{
				seq.insert(seq.end(), line_begin, line_begin + fasta.gcount() - 1);
				break;
			}
			seq.insert(seq.end(), line_begin, line_end);
		}
		seq.shrink_to_fit();
		return new FastaRecord(std::move(name), std::move(seq));
	}
}

/*
 * Data la descrizione del cromosoma se ne estrae il nome
 */
void ParseFasta::ExtractNameFromDescription(std::string& name)
{
	std::size_t found = name.find(" ");
	if(found !=  std::string::npos)
		name = name.substr (1,found - 1);
	else
		name = name.substr (1,name.size() - 1);

	if(name.find("chr") == std::string::npos)
		name = "chr" + name;
}


/*
 * Questo metodo popola il vettore chromosome il quale è un vettore di FastaRecord(struttura di
 * due stringhe una contenente il nome, l'altra la sequenza
 */
void ParseFasta::read_fasta(const std::string& filename_fasta)
{
	std::ifstream fasta(filename_fasta);

	size_t num_records = count_records(fasta, fasta_delim);

	//std::cerr <<"num records"<<num_records << std::endl;
	chromosome.reserve(num_records);

	for (; num_records > 0; --num_records)
	{
		FastaRecord* currentChromosome = read_fasta_record(fasta);
		number_of_genome_bases += currentChromosome->sequence.size();
		ExtractNameFromDescription(currentChromosome->name);
		chromosome.emplace_back(currentChromosome);
	}
}

/*
 *
 */
FastaRecord* ParseFasta::GetSequence(size_t i)
{
	return chromosome[i];
}

/**
 *
 */
size_t ParseFasta::Size()
{
	return chromosome.size();
}

/**
 *
 */
FastaRecord* ParseFasta::find_seq_by_name(const std::string current_chr_name, unsigned int& start_index_chr)
{
	//For chromosome 1 the start index is 1
	unsigned int start_index_true = 1;

	//For each chromosome in which the reference genome is divided
	for (size_t i = 0; i < chromosome.size(); i++)
	{
		//Get the sequence name
		const std::string seq_name = GetSequence(i)->name;

		//Find the chromosome containing the current analyze SNP
		if (current_chr_name.compare(seq_name) == 0)
		{
			start_index_chr = start_index_true;

			//return name,sequence and size of chromosome
			return GetSequence(i);
		}
		else
			//For chromosome 2 start index is size of chromosome 1 etc
			start_index_true += GetSequence(i)->sequence.size();
	}

	//If the current snp is contained in a chromosome not contained in the reference genome
	start_index_chr = 0;
	return nullptr;
}

uint64_t ParseFasta::GetNumOfGenomeBases()
{
	return number_of_genome_bases;
}

bool IsChromosomePresent(std::string chromosome_name)
{
	return false;
}
