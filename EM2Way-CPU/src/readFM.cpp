/*
 * readFM.cpp
 *
 *  Created on: Jan 31, 2019
 *      Author: hn915
 */


#include "readFM.h"

void read_meta(FILE * FM_meta_fp, 
                uint64_t * fmt_len, 
                bool * c32, 
                uint32_t * bucket_bwt_len,
		        uint64_t * endCharPos, 
                uint32_t * bucket_pad_size) {

    uint32_t BP_bit;
    uint32_t BP_range;
    uint32_t bucket_size;

    readFile(FM_meta_fp, endCharPos, sizeof(uint64_t));

    readFile(FM_meta_fp, fmt_len, sizeof(uint64_t));

    readFile(FM_meta_fp, &BP_bit, sizeof(uint32_t));
    //check if the meta-data is the same as the index
    if (FM_BP_BIT != BP_bit) {
        fprintf(stderr, "error: FM_BP_RANGE is different from the index!\n");
        exit(1);
    }

    readFile(FM_meta_fp, &BP_range, sizeof(uint32_t));
    //check if the meta-data is the same as the index
    if (FM_BP_RANGE != BP_range) {
        fprintf(stderr, "error: FM_BP_RANGE is different from the index!\n");
        exit(1);
    }

    readFile(FM_meta_fp, &bucket_size, sizeof(uint32_t));
    //check if the meta-data is the same as the index
    if (BUCKET_SIZE != bucket_size) {
        fprintf(stderr, "error: BUCKET_SIZE is different from the index!\n");
        exit(1);
    }

    readFile(FM_meta_fp, c32, sizeof(bool));
    readFile(FM_meta_fp, bucket_bwt_len, sizeof(uint32_t));

    readFile(FM_meta_fp, bucket_pad_size, sizeof(uint32_t));

}
