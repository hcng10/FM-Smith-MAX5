#ifndef WRITESAM_H
#define WRITESAM_H

#include <stdint.h>
#include <omp.h>

#include "reads.h"
#include "format.h"
#include "fm-sa-def.h"

#define POS_BIT 10
#define QS_BIT 3
#define MATCH_BIT 5
#define ED_DIST_BIT 2

#define ALIGN_REPORT_NUM 16

void convertSAM(FILE *fp_sam, std::vector<read2Bit_t> &reads, uint32_t *sai, char *buffer, std::vector<chr_t> &chrs);

#endif