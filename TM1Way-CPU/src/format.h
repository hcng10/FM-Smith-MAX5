#ifndef FM_SA_FORMAT_H
#define FM_SA_FORMAT_H


#include <cstdlib>
#include <string>
#include <vector>
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

void readNinfo(FILE * fp, uint64_t N_cluster, std::vector<nchar_cluster_t> &nchar_clusters);
void readChrName(FILE * fp, uint16_t chrs_num, std::vector<chr_t> &chrs);

#endif //FM_SA_H