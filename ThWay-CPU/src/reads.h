/*
 * reads.h
 *
 *  Created on: Feb 3, 2019
 *      Author: hn915
 */

#ifndef READS_H_
#define READS_H_

#include <vector>
#include <queue>
#include <omp.h>
#include <stdio.h>
#include "file_op.h"
#include <sys/time.h>

struct seedhit_t {
	uint8_t seg;
	uint32_t low;
	uint32_t high;
    uint32_t range;
    bool is_rev_cmpt;
};

struct SARamgeCmp{

    bool operator()(const seedhit_t &s_a, 
                    const seedhit_t &s_b){

        return s_a.range > s_b.range;   
    }
};


struct read2Bit_t {
    char seq[MAX_READ_LEN+1];
    uint32_t seq_len;

    // symbol in binary
    uint8_t pck_sym[CEIL(MAX_READ_LEN, FM_BP_RANGE)];
    uint8_t pck_sub_sym[MISMATCH_DIV][CEIL(MAX_READ_LEN, FM_BP_RANGE)/MISMATCH_DIV + 1];
    uint8_t seq_sub_len[MISMATCH_DIV];
    uint8_t seq_sub_pos[MISMATCH_DIV];

    std::string at_line;
    std::string q_score;

    //char plus_data;
    bool is_align;
    uint32_t n_seedhits;

    priority_queue<seedhit_t, vector<seedhit_t>, SARamgeCmp> seedhit;
    //std::vector<seedhit_t> seedhit;
};


struct read3Bit_t {
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
};

void loadReads(FILE *fp, std::vector<read2Bit_t> &reads2Bit, 
		char *buffer, uint64_t size_r, uint64_t *bytes, bool r_ctrl);

uint64_t writeReads(FILE *fp, std::vector<read2Bit_t> &reads, char *buffer);
uint64_t writeReadsN(FILE *fp, std::vector<read3Bit_t> &reads, char *buffer);


#endif /* READS_H_ */
