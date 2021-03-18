//
// Created by hcng on 8/9/18.
//

#ifndef FM_SA_H
#define FM_SA_H

#define HARD_CODE_READ_LEN 101
// tmp buf size
#define LINE_SIZE 512
#define MAX_RD_CNT 1000000
#define BUFF_SIZE 800000000

#define BP_BIT 3
#define QS_BIT 6

#define PADSIZE 72
#define RSLT_PADSIZE 24
#define RSLT_TAIL_ALIGN 1

#define REF_DEF_LEN 160
#define READ_DEF_LEN 104

#define SEED_LEN 10
#define TRIM_READ_LEN 16

#define MAX_EXTD 16


#define N_KRNL 1
#define LATENCY 25

#define IS_SIM 1


//#define MAX_READ_LEN 2000

//800000000

// round up division
#define CEIL(a, b) (((a)+(b)-1)/(b))

#define ALIGN_REPORT_NUM 16
#define EXTEND_PER_SEED 15

#endif //FM_SA_H
