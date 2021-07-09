#ifndef GENO_H
#define GENO_H

#include <string>
#include <iostream>
#include <cstdint>
#include "loadSmartDictionary.h"
#include "utility.h"
#include <math.h>
#include <unordered_map>
#include <unordered_set>
#include "fq.hpp"
#include "genolight.h"

#define MAX_COV ((1 << 6) - 1)

enum
{
	GTYPE_NONE, GTYPE_REF, GTYPE_ALT, GTYPE_HET
};

struct chrlens
{
		std::string name;
		size_t len;

		chrlens(std::string name, size_t len) :
		name
		{name}, len
		{len}
		{
		}
};

struct call
{
		int genotype;
		double confidence;
};

#define CALL(g, c) ((struct call){.genotype = (g), .confidence = (c)})

void genotype(param_creation_geno* geno_param);

inline struct call choose_best_genotype(const int ref_cnt, const int alt_cnt, const uint8_t ref_freq_enc, const uint8_t alt_freq_enc);

param_creation_geno geno_param_init(int argc, char *argv[]);

#endif
