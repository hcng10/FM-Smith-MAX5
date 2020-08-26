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


struct read2Bit_t {
    char seq[MAX_READ_LEN+1];
    uint32_t seq_len;

    // symbol in binary
    uint8_t pck_sym[CEIL(MAX_READ_LEN, FM_BP_RANGE)];

    std::string at_line;
    std::string q_score;

    //char plus_data;
    uint32_t low_bw;
    uint32_t high_bw;

    uint32_t low_fw;
    uint32_t high_fw;


    bool isaligned_bw;
    bool isaligned_fw;
    bool has_N;
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

void loadReads(FILE *fp, std::vector<read2Bit_t> &reads2Bit, std::vector<read3Bit_t> &reads3Bit,
		char *buffer, uint64_t size_r, uint64_t *bytes, bool r_ctrl);

uint64_t writeReads(FILE *fp, std::vector<read2Bit_t> &reads, char *buffer);
uint64_t writeReadsN(FILE *fp, std::vector<read3Bit_t> &reads, char *buffer);


#endif /* READS_H_ */
