/*
 * align.h
 *
 *  Created on: Feb 6, 2019
 *      Author: hn915
 */

#ifndef ALIGN_H_
#define ALIGN_H_


#include <vector>
#include <mutex>   
#include <sys/time.h>
#include <MaxSLiCInterface.h>


#include "readFM.h"
#include "reads.h"
#include "fm-sa-def.h"





// kernel input
struct in_t {
  uint32_t id;
  uint8_t pck_sym[CEIL(MAX_READ_LEN, FM_BP_RANGE)];
  uint8_t len;
  uint8_t is_pad;
};

// kernel output
struct out_t {
  uint32_t id;
  uint32_t low_bw;
  uint32_t high_bw;
  uint32_t low_fw;
  uint32_t high_fw;
  uint32_t isaligned_bw;
  uint32_t isaligned_fw;
  uint32_t pad;
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
	double parse_time;
	uint64_t aligned_cnt;
};



/*void writeIndexBothWay(index32_t * index32, index32_t * index32_rev,
						uint64_t index_bytes,
						uint32_t n_buckets,
						fpga_t fpga_var);*/

max_run_t *writeIndex(index32_t * index, uint64_t index_bytes, fpga_t fpga_var);

//align reads
void align(std::vector<read2Bit_t> &reads,
		uint64_t index_bytes,
		uint32_t high_init,
		uint32_t bucket_bwt_len,
		uint32_t end_char_pos,
		fpga_t fpga_var,
		bmk_t * bmk);



#endif /* ALIGN_H_ */
