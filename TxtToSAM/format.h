#ifndef FM_SA_FORMAT_H
#define FM_SA_FORMAT_H


#include <cstdlib>
#include <string>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "fm-sa-def.h"

#define NAME_LEN 1000

using namespace std;

struct nchar_cluster_t{
    uint64_t ref_cnt;
    uint64_t fmt_cnt;
    uint32_t un_cnt;
    uint64_t cum_un_cnt;
};

struct chr_t{
    string name;
    uint64_t begin;
    uint64_t end;
};

struct hit_om_t {
	uint8_t mis_pos;
	uint32_t low;
	uint32_t high;
	char mis_sym;
	char qs;
    bool is_aligned_bck;
    bool is_rev_cmpt;
};

struct hit_tm_t {
	uint8_t mis_pos0;
	uint8_t mis_pos1;
	uint32_t low;
	uint32_t high;
	char mis_sym0;
	char mis_sym1;
	char qs0;
	char qs1;
    bool is_rev_cmpt;
};


struct read2Bit_t {
    uint32_t seq_len;
    std::string seq;

    std::string at_line;
    std::string q_score;

    uint8_t n_hits;
    
    std::vector<hit_om_t> hits_om;
    std::vector<hit_tm_t> hits_tm;
};


void readNinfo(FILE * fp, uint64_t N_cluster, std::vector<nchar_cluster_t> &nchar_clusters);
void readChrName(FILE * fp, uint16_t chrs_num, std::vector<chr_t> &chrs);

void loadtxt(FILE *fp, 
                std::vector<read2Bit_t> &reads2Bit, 
                bool is_om,
                char *buffer, 
                uint64_t size_r, 
                uint64_t *bytes, 
                bool r_ctrl);

#endif //FM_SA_H