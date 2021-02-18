#include "format.h"

void readNinfo(FILE * fp, uint64_t N_cluster, std::vector<nchar_cluster_t> &nchar_clusters){

    nchar_cluster_t nchar_cluster;

    for (uint64_t i = 0; i < N_cluster; i++){
        fread(&(nchar_cluster.ref_cnt), sizeof(char), sizeof(uint64_t), fp);
        fread(&(nchar_cluster.fmt_cnt), sizeof(char), sizeof(uint64_t), fp);
        fread(&(nchar_cluster.un_cnt), sizeof(char), sizeof(uint32_t), fp);
        fread(&(nchar_cluster.cum_un_cnt), sizeof(char), sizeof(uint64_t), fp);

        nchar_clusters.push_back(nchar_cluster);
    }
}

void readChrName(FILE * fp, uint16_t chrs_num, std::vector<chr_t> &chrs){

    chr_t chr;
    uint64_t length;
    char name[NAME_LEN];

    for (uint16_t i = 0; i < chrs_num; i++){
        fread (&(chr.begin), sizeof(uint64_t), 1, fp);
        fread(&(chr.end), sizeof(uint64_t), 1, fp);
        fread(&length, sizeof(uint64_t), 1, fp);
        fread(&name, sizeof(char), length, fp);

        name[length] = '\0';

        //chr.name = name;
        chr.name = string(name, length); 

        chrs.push_back(chr);
    }
}