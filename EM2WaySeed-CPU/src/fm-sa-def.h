/*
 * bowtie2-lite-def.h
 *
 *  Created on: Jan 30, 2019
 *      Author: hn915
 */

#ifndef FM_SA_DEF_H
#define FM_SA_DEF_H


#define BUCKET_SIZE 64
#define FM_BP_BIT 2
#define FM_BP_RANGE 4

#define FM_BPN_BIT 3
#define FM_BPN_RANGE 5

#define MAX_READ_LEN 168
#define BUFF_SIZE 4000000000
//160000000

#define N_DFE 1
#define N_KRNL 3
#define BURST_BYTES 64
#define LATENCY 384


// round up division
#define CEIL(a, b) (((a)+(b)-1)/(b))

//related to seeding 
#define ACTUAL_READ_LEN 101
#define SEED_INTERVAL 8
#define SEED_LEN 20 //it is fixed in Bowtie2

#define SEED_NUM (ACTUAL_READ_LEN-(SEED_LEN-SEED_INTERVAL))/SEED_INTERVAL
#define SEED_NUM_CEIL SEED_NUM+1
#define IS_REV_CMPT_BIT 2

#define HIT_PAD 9

#define EXTEND_PER_SEED 16

#define IS_SIM 0


#endif /* FM_SA_DEF_H */

#include "EM2WaySeed.h"
