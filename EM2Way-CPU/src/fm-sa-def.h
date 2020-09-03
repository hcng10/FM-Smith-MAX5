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
#define BUFF_SIZE 8000000000
//160000000

#define N_DFE 1
#define N_KRNL 3
#define BURST_BYTES 64
#define LATENCY 384


// round up division
#define CEIL(a, b) (((a)+(b)-1)/(b))



#endif /* FM_SA_DEF_H */

#include "EM2Way.h"
