#ifndef VARGENO_H
#define VARGENO_H

#include <stdint.h>
#include <atomic>

/////////////////////////////
//#define REF_LITE     0
//#define PCOMPACT     0

#define ERR_RATE       0.01
#define AVG_COV        7.1
#define MAX_MATES_DIST 2000
/////////////////////////////

#define BASE_A 0x0
#define BASE_C 0x1
#define BASE_G 0x2
#define BASE_T 0x3
#define BASE_N 0x4
#define BASE_X 0x7

/* Note:
 * Frequencies (between 0 and 1) are encoded as 8-bit uints as: freq*0xff.
 * They are decoded as: freq_enc/255.0f
 */

struct snp_kmer_entry
{
		uint64_t kmer_lo40 :40;
		snp_info snp;
		uint32_t pos;
		uint8_t ambig_flag;
}__attribute__((packed));

#if !PCOMPACT
struct packed_pileup_entry
{
		unsigned ref :2;
		unsigned alt :2;
		//unsigned ref_cnt :6;
		//unsigned alt_cnt :6;
		std::atomic<int_least8_t> ref_cnt;
		std::atomic<int_least8_t> alt_cnt;
		uint8_t ref_freq;
		uint8_t alt_freq;
}__attribute__((packed));
#endif

#endif

