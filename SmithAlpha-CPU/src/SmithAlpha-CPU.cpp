#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <stdint.h>
#include <bitset>
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <math.h> 
#include <thread>
#include <omp.h>

#include "fm-sa-def.h"

#include "file_op.h"
#include "format.hpp"
#include "align.hpp"
#include "aligned_rslt.hpp"

#include "Maxfiles.h"

using namespace std;

void loadseeds(FILE *seed_fp, 
            vector<Read_t> & read_vec,
            vector<InToRead_t> & intoread_vec,
            char * ref_fmt, 
            uint64_t ref_fmt_len, 
            uint64_t & total_extend_cnt){


    rdSeed(seed_fp, read_vec); 
    packRdRef(read_vec, intoread_vec, ref_fmt, ref_fmt_len, total_extend_cnt);

    return;
}

int main(int argc, char *argv[]) {

    //struct timeval  tv1, tv2;

#if IS_SIM == 0
    // program usage
    if (argc != 4) {
        printf("usage: %s <FASTA file> <Seed file> <Output extend file>\n", argv[0]);
        exit(1);
    }

    string s_in = string(argv[1]);
    string seed_in = string(argv[2]);
    //string extend_out = string(argv[3]);

#elif IS_SIM == 1

    string s_in = "sample_index/sample_922char_wN.fa";
    string seed_in = "sample_index/em2w.txt";
    //string extend_out = "sample_index/extend_out.txt";
#endif


    //string r1 = "sm.fq";
    string ra1 = "sm.txt";


    char * out_buff = new char [BUFF_SIZE];
    char * aligned_buff = new char [BUFF_SIZE];

    // the input refernce fasta file
    FILE * fasta_fp = NULL;

    // temp storage for the reference (ref), in char 
    char * fmt = NULL;
    // convert the nucleotide to 4 bits, pack them
    uint8_t * pck_fmt = NULL;


    // vector containing chromosome name
    uint16_t chrs_cnt = 0;
    std::vector<chr_t> chrs;
    
    

    // 1. construct reference sequence in binary format
    printf("Getting reference sequence ... \n"); fflush(stdout);

    openFile(&fasta_fp, s_in, "r");
    uint64_t fa_len = fileSizeBytes(fasta_fp);//uint64_t 

    uint64_t fmt_len = 0;
    //uint64_t ref_len = 0;
    
       // new char array to hold all the nucleotide, 1 byte 1 nucleotide
    fmt = new char [fa_len];
    if (!fmt){
        fprintf(stderr, "error: unable to allocate memory for reference sequence!\n");
        exit(1);
    }

    faToFmt(fasta_fp, fmt, fmt_len, chrs);

    printf("FINISH ---> Getting reference\n\n");fflush(stdout);


    

    // 2. program FPGA
    printf("\nInitializing FPGA ... \n"); fflush(stdout);
    fpga_t fpga_var;
	
    fpga_var.maxfile = SmithAlpha_init();
	fpga_var.group = max_load_group(fpga_var.maxfile, MAXOS_EXCLUSIVE, "*", 1);
	fpga_var.engine = max_lock_any(fpga_var.group);

    printf("OK\n"); fflush(stdout);

    

    
    // 3. read in the seed info
    FILE * seed_fp = NULL;
    openFile(&seed_fp, seed_in, "r");

    std::vector<Read_t> read_vec1;
    std::vector<Read_t> read_vec2;

    std::vector<InToRead_t> intoread_vec1;
    std::vector<InToRead_t> intoread_vec2;

    // 3. file to write aligned result
    FILE *aligned_fp = NULL;
    openFile(&aligned_fp, ra1, "w+");
    

    uint64_t total_extend_cnt1 = 0;
    uint64_t total_extend_cnt2 = 0;

    // read the READs
    printf("Start reading one batch of read ...\n");
    
    loadseeds(seed_fp,
                read_vec1,
                intoread_vec1,
                fmt,
                fmt_len,
                total_extend_cnt1);


    //rdSeed(seed_fp, read_vec1);
    //packRdRef(read_vec1, fmt, fmt_len, total_extend_cnt1);

    //cout<<"Final dai jei "<<total_extend_cnt1<<"\n";

    bmk_t bmk1, bmk2;
	bmk1.in_size = 0;
	bmk1.out_size = 0;
	bmk1.process_time = 0;
	bmk1.aligned_cnt = 0;

	bmk2.in_size = 0;
	bmk2.out_size = 0;
	bmk2.process_time = 0;
	bmk2.aligned_cnt = 0;

    uint32_t cnt = 0;

    std::thread loadseeds_thread;
    std::thread writeAligned_thread;

    // 1 thread to read, 1 thread to align on FPGA
    for (char thread_cnt = 0; ; thread_cnt++) {

        // load to reads2, process reads1
        if (!(thread_cnt % 2)) {

            printf("\nLoad Read (A) Again\n"); fflush(stdout);
            loadseeds_thread = std::thread(loadseeds,
                                        seed_fp,
                                        std::ref(read_vec2),
                                        std::ref(intoread_vec2),
                                        fmt,
                                        fmt_len,
                                        std::ref(total_extend_cnt2));
            /*loadseeds(seed_fp,
                    std::ref(read_vec2),
                    std::ref(intoread_vec2),
                    fmt,
                    fmt_len,
                    std::ref(total_extend_cnt2));*/
            
            if (read_vec1.size() > 0) {
                cnt += read_vec1.size();


                align(read_vec1, intoread_vec1, total_extend_cnt1, fpga_var, &bmk1);
                total_extend_cnt1 = 0;


                writeAligned_thread = std::thread(writePosDirectly, 
                                aligned_fp, 
                                std::ref(read_vec1),
                                std::ref(intoread_vec1),
                                std::ref(chrs),
                                aligned_buff);
                
                writeAligned_thread.join();


            }else{
                total_extend_cnt1 = 0;
                break;
            }
        }else{

            printf("\nLoad Read (B) Again\n"); fflush(stdout);
            loadseeds_thread = std::thread(loadseeds,
                                        seed_fp,
                                        std::ref(read_vec1),
                                        std::ref(intoread_vec1),
                                        fmt,
                                        fmt_len,
                                        std::ref(total_extend_cnt1));
            
            if (read_vec2.size() > 0) {
                cnt += read_vec2.size();

                align(read_vec2, intoread_vec2, total_extend_cnt2, fpga_var, &bmk2);
                total_extend_cnt2 = 0;


            }else{
                total_extend_cnt2 = 0;
                break;
            }

        }
        loadseeds_thread.join();    
    
    }
    loadseeds_thread.join();

    printf("processed %u reads\n", cnt);

    delete [] out_buff;
    delete [] aligned_buff;
    delete [] fmt;
    
}
