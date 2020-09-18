/*
 * reads.h
 *
 *  Created on: Feb 3, 2019
 *      Author: hn915
 */

#ifndef READS_H_
#define READS_H_

#include <vector>
#include <omp.h>
#include <stdio.h>
#include "file_op.h"
#include <sys/time.h>

// hit
/*struct hit_t {
	uint8_t mis_pos[N_HITS];
	uint32_t low[N_HITS];
	uint32_t high[N_HITS];
	char mis_sym[N_HITS];
	uint8_t qs[N_HITS];

};*/

struct hit_t {
	uint8_t mis_pos_sorted[N_HITS];
	uint32_t low_sorted[N_HITS];
	uint32_t high_sorted[N_HITS];
	char mis_sym_sorted[N_HITS];
	char qs_sorted[N_HITS];
    uint32_t is_aligned_bck;
	uint8_t id;
	uint8_t n_hits;
    bool is_rev_cmpt;
	//uint8_t pad;
};


struct read2Bit_t {
    char seq[MAX_READ_LEN+1];
    uint32_t seq_len;

    // symbol in binary
    uint8_t pck_sym[CEIL(MAX_READ_LEN, FM_BP_RANGE)];
    uint8_t pck_sym_qs[MAX_READ_LEN];//2bit sym, 6bit qs

    std::string at_line;
    std::string q_score;

    //char plus_data;
    //uint32_t low;
    //uint32_t high;

    bool is_align;
    uint8_t n_hits;
    //bool has_N;

    std::vector<hit_t> hits;
};


/*struct read3Bit_t {
    char seq[MAX_READ_LEN+1];
    uint32_t seq_len;

    // symbol in binary
    uint8_t pck_sym[CEIL(MAX_READ_LEN, FM_BPN_RANGE)];

    std::string at_line;
    std::string q_score;

    //char plus_data;
    uint32_t low;
    uint32_t high;

    bool is_align;
    bool has_N;
};*/

void loadReads(FILE *fp, std::vector<read2Bit_t> &reads2Bit, //std::vector<read3Bit_t> &reads3Bit,
		char *buffer, uint64_t size_r, uint64_t *bytes, bool r_ctrl);

uint64_t writeReads(FILE *fp, std::vector<read2Bit_t> &reads, char *buffer);
//void writeReadsN(FILE *fp, std::vector<read3Bit_t> &reads, char *buffer);


#endif /* READS_H_ */
