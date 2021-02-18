/*
 * align.h
 *
 *  Created on: Feb 6, 2019
 *      Author: hn915
 */

#ifndef ALIGN_H_
#define ALIGN_H_

#include <vector>
#include <sys/time.h>
#include <thread>
#include <MaxSLiCInterface.h>


#include "readFM.h"
#include "reads.h"
#include "format.h"
#include "fm-sa-def.h"
#include "Maxfiles.h"




// kernel input
struct in_t {
  uint32_t id;
  //uint8_t pck_sym[CEIL(MAX_READ_LEN, FM_BP_RANGE)];
  uint8_t pck_sym_qs[MAX_READ_LEN];
  uint8_t len;
  uint8_t is_pad;
};

// kernel output
/*struct out_t {
  uint32_t id;
  uint32_t low;
  uint32_t high;
  uint32_t pad;
};*/

// kernel output

struct out_t {
	uint8_t mis_pos1_sorted[N_HITS];
	uint8_t mis_pos2_sorted[N_HITS];
	uint32_t low_sorted[N_HITS];
	uint32_t high_sorted[N_HITS];
	uint8_t qs_mis_sym1_sorted[N_HITS];
	uint8_t qs_mis_sym2_sorted[N_HITS];
	uint32_t is_aligned_bck;
	//uint32_t is_rev_cmpt; // is a bitwise array
	uint32_t id;
	uint8_t n_hits;
	uint8_t pad[(HIT_PAD)/(sizeof(uint8_t)*8)];
};



struct fpga_t{
	max_file_t * maxfile;
	max_engine_t * engine;
	max_group_t * group;
};

struct bmk_t{
	uint64_t in_size;
	uint64_t out_size;
	double process_time;
	uint64_t aligned_cnt;
};



max_run_t *writeIndex(index32_t * index, uint64_t index_bytes, fpga_t fpga_var);
void writeIndexBothWay(index32_t * index, index32_t * index_rev, uint64_t index_bytes, fpga_t fpga_var);

//align reads
void align(std::vector<read2Bit_t> &reads,
			//uint32_t * cnt32,
			//uint32_t * cnt32_rev,
			uint64_t index_bytes,
			//uint64_t index_bytes_rev,
			uint32_t high_init,
			uint32_t high_init_rev,
			uint32_t bucket_bwt_len,
			uint32_t bucket_bwt_len_rev,
			uint32_t end_char_pos,
			uint32_t end_char_pos_rev,
			fpga_t fpga_var,
			bmk_t * bmk);

/*void align(std::vector<read2Bit_t> &reads,
			//uint32_t * cnt32,
			//uint32_t * cnt32_rev,
			uint64_t index_bytes,
			//uint64_t index_bytes_rev,
			uint32_t high_init,
			//uint32_t high_init_rev,
			uint32_t bucket_bwt_len,
			//uint32_t bucket_bwt_len_rev,
			uint32_t end_char_pos,
			//uint32_t end_char_pos_rev,
			fpga_t fpga_var,
			bmk_t * bmk);*/


#endif /* ALIGN_H_ */
