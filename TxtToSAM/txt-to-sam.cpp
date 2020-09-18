#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <string>
#include <cstring>
#include <bitset>
#include <thread>
#include <functional> 

#include "format.h"
#include "readFM.h"
#include "file_op.h"
#include "writeSAM.hpp"


using namespace std;

int main(int argc, char *argv[]) {

    // stands for reconfigurable alignment
    string ext = "ral";

    // file descriptor for the input files
    FILE * N_fp = NULL;
    FILE * FM_fp = NULL;
    FILE * FM_meta_fp = NULL;

    FILE * rev_N_fp = NULL;
    FILE * rev_FM_fp = NULL;
    FILE * rev_FM_meta_fp = NULL;

    FILE * SA_fp = NULL;
    FILE * rev_SA_fp = NULL;

    FILE * chr_fp = NULL;

    FILE * txt_fp = NULL;
    FILE * sam_fp = NULL;

    uint64_t f_size;

    // fmt_len is the length of the reference without the N nucleotide
    uint64_t fmt_len;
    bool is32bit;

    uint32_t bucket_bwt_len;
    uint32_t bucket_pad_size;

    uint64_t end_char_pos;

    // suffix array
    uint32_t * sai = NULL;
    uint32_t * sai_rev = NULL;

    // number of chromosomes
    vector<chr_t> chrs;
    uint16_t chrs_num;
    
    // N char info
    vector<nchar_cluster_t> nchar_clusters;
    uint64_t N_cluster;

    bool is_om = false;

    // program usage
    if (argc != 5) {
        printf("usage: %s <index basename> <text file> <SAM file> <is_OM>\n", argv[0]);
        exit(1);
    }

    string s1 = string(argv[1]) + ".1." + ext;
    string s2 = string(argv[1]) + ".2." + ext;
    string s3 = string(argv[1]) + ".3." + ext;

    string s4 = string(argv[1]) + ".4." + ext;
    string s5 = string(argv[1]) + ".5." + ext;
    string s6 = string(argv[1]) + ".6." + ext;

    string s7 = string(argv[1]) + ".7." + ext;
    string s8 = string(argv[1]) + ".8." + ext;

    string s9 = string(argv[1]) + ".9." + ext;
    string s10 = string(argv[1]) + ".10." + ext;

    cout<<"Read file name: "<<argv[1]<<"\n";

    is_om =  strcmp(argv[4], "true") == 0? true: false;
    
    // load index
    printf("Reading meta data ... "); fflush(stdout);

    openFile(&FM_fp, s1, "r");
    openFile(&FM_meta_fp, s2, "r");


    f_size = fileSizeBytes(FM_fp);
    fclose(FM_fp);

    read_meta(FM_meta_fp, 
                &fmt_len, 
                &is32bit, 
                &bucket_bwt_len,
                &end_char_pos, 
                &bucket_pad_size,
                &N_cluster,
                &chrs_num);


    fclose(FM_meta_fp);
    printf("FINISH ---> Reading meta data\n\n");fflush(stdout);


    // read SA
    printf("Reading SA ... \n");fflush(stdout);

    openFile(&SA_fp, s7, "r");
    
    sai = new uint32_t[fmt_len];

    for (uint32_t i = 0; i < fmt_len; i++){
        size_t sizeread = fread(&sai[i], 1, sizeof(uint32_t), SA_fp);
        if (sizeread != sizeof(uint32_t)) {
            fprintf(stderr, "error: unable to read SA file!\n");
            exit(1);
        }
    }

    fclose(SA_fp);
    printf("FINISH ---> Reading SA\n\n");


    // read SA
    if (is_om == true){
        printf("Reading Reversed SA ... \n");fflush(stdout);

        openFile(&rev_SA_fp, s8, "r");
        
        sai_rev = new uint32_t[fmt_len];

        for (uint32_t i = 0; i < fmt_len; i++){
            size_t sizeread = fread(&sai_rev[i], 1, sizeof(uint32_t), rev_SA_fp);
            if (sizeread != sizeof(uint32_t)) {
                fprintf(stderr, "error: unable to read SA file!\n");
                exit(1);
            }
        }

        fclose(rev_SA_fp);
        printf("FINISH ---> Reading SA\n\n");
    }


    // read the name of chromosome
    printf("Reading the name of chromosome ... \n");fflush(stdout);

    openFile(&chr_fp, s9, "r");
    readChrName(chr_fp, chrs_num, chrs);

    fclose(chr_fp);
    printf("FINISH ---> Reading name of chromosome\n\n");


    //read text file
    openFile(&txt_fp, argv[2], "r");
    openFile(&sam_fp, argv[3], "w+");


    // allocate I/O buffers
    char * in_buff = new char [BUFF_SIZE + 512];
    char * sam_buff = new char [BUFF_SIZE];

    std::vector<read2Bit_t> reads1, reads2;
    f_size = fileSizeBytes(txt_fp);

    // determine how much data to read from the file
    uint64_t bytes_r = 0;
    uint64_t size_r =  bytes_r + BUFF_SIZE <= f_size ? BUFF_SIZE : f_size - bytes_r;

    std::thread loadread_thread;
    cout<<f_size<<" "<<bytes_r<<" "<<size_r<<"\n";
    printf("Load first batch of result ... \n");fflush(stdout);
    loadtxt(txt_fp, reads1, is_om, in_buff, size_r, &bytes_r, true);
    cout<<bytes_r<<" "<<size_r<<"\n";

    for (char thread_cnt = 0; ; thread_cnt++) {
     	bool r_ctrl = bytes_r < f_size ? true : false;
    	size_r =  bytes_r + BUFF_SIZE <= f_size ? BUFF_SIZE : f_size - bytes_r;   

    cout<<bytes_r<<" "<<size_r<<"\n";
        // load to reads2, process reads1
        if (!(thread_cnt % 2)) {
            printf("\nLoad Result (A) Again\n"); fflush(stdout);

            loadread_thread = std::thread(loadtxt,
                                            txt_fp,
                                            std::ref(reads2),
                                            is_om,
                                            in_buff,
                                            size_r, 
                                            &bytes_r, 
                                            true);   

            if (reads1.size() > 0) {
                convertSAM(sam_fp, sam_buff, reads1, sai, sai_rev, chrs, is_om);

                printf("FINISH: One round conversion: count %d\n", reads1.size());

                printf("...................One pass finished.......................................\n"); fflush(stdout);

            }else{
                break;
            }
        }else{
            printf("\nLoad Result (B) Again\n"); fflush(stdout);

            loadread_thread = std::thread(loadtxt,
                                            txt_fp,
                                            std::ref(reads1),
                                            is_om,
                                            in_buff,
                                            size_r, 
                                            &bytes_r, 
                                            true);   
            if (reads2.size() > 0) {
                convertSAM(sam_fp, sam_buff, reads2, sai, sai_rev, chrs, is_om);

                printf("FINISH: One round conversion: count %d\n", reads2.size());

                printf("...................One pass finished.......................................\n"); fflush(stdout);
            }else{
                break;
            }
        }
        loadread_thread.join();
    }
    loadread_thread.join();

    fclose(txt_fp);
    fclose(sam_fp);

    delete [] sai;
    if (sai_rev != NULL) delete [] sai_rev;
    delete [] in_buff;
    delete [] sam_buff;

    //convertSAM(sam_fp, sam_buff, reads1, sai, sai_rev, chrs, is_om);
}