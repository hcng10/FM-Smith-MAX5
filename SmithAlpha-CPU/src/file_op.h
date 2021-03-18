//
// Created by hcng on 6/27/18.
//

#ifndef FM_SA_FILE_OP_H
#define FM_SA_FILE_OP_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <bitset>
#include <vector>
#include <queue>

#include "fm-sa-def.h"


using namespace std;

// open file
void openFile(FILE **fp, string f_name, const char *mode);

/**
    Uses fstat.st_size to get the file size in byte

    @param  *fp     input fd for the fasta file

    @return file size in byte
*/
uint64_t fileSizeBytes(FILE *fp);

void writeFile(FILE *fp, void *data, uint64_t n_bytes);

/**
    Writes info for a N nucleotide cluster with the following
    (A) the starting point of N nucleotides in the original reference (64bit)
    (B) the point where N nucleotide got cut off and replaced 
    with actual nucleotie in the reconstructed reference (64bit)
    (C) Number of N nucleotides in this cluster (32bit)

    @param  *fp         input fd for the fasta file
    @param  ref_cnt     the starting point of N nucleotide in original reference
    @param  fmt_cnt     the point where N nucleotide got cut off in the reconstructed ref
    @param  un_cnt      number of N nucleotide in this cluster
    @param  cum_un_cnt  cumulative number of N char up to this point 

*/
void writeNinfo(FILE * fp, uint64_t ref_cnt, uint64_t fmt_cnt, uint32_t un_cnt, uint64_t cum_un_cnt);

#endif //FM_SA_H
