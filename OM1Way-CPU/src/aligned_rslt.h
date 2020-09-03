#ifndef ALIGNED_RSLT_H_
#define ALIGNED_RSLT_H_

#include <vector>
#include <omp.h>
#include <stdio.h>
#include <sys/time.h>


#include "file_op.h"
#include "reads.h"

#define POS_BIT 10
#define QS_BIT 3
#define MATCH_BIT 5
#define ED_DIST_BIT 2

void writeAligned(FILE *fp, std::vector<read2Bit_t> &reads, char *buffer);

#endif