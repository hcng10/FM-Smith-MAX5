/*
 * readFM.h
 *
 *  Created on: Jan 30, 2019
 *      Author: hn915
 */

#ifndef READFM_H_
#define READFM_H_

#include "fm-sa-def.h"
#include "file_op.h"
#include <stdint.h>

#define BWT32_LEN_wPADDED ((BUCKET_SIZE * 8 - 32 * FM_BP_RANGE) / 8)


struct index32_t {
    uint32_t count[FM_BP_RANGE];
    uint8_t bwt[BWT32_LEN_wPADDED];

};

void read_meta(FILE * FM_meta_fp, 
                uint64_t * fmt_len, 
                bool * c32, 
                uint32_t * bucket_bwt_len,
		        uint64_t * endCharPos, 
                uint32_t * bucket_pad_size,
                uint64_t * N_cluster,
                uint16_t * chrs_num);

#endif /* READFM_H_ */
