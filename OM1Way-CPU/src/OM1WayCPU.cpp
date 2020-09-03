#include <stdio.h>
#include <stdint.h>
#include <string>
#include <cstring>
#include <iostream>
#include <bitset>
#include <thread>

#include "Maxfiles.h"
#include "reads.h"
#include "readFM.h"
#include "file_op.h"
#include "align.h"
#include "aligned_rslt.h"
#include "fm-sa-def.h"

#define IS_SIM 1

using namespace std;

template<typename T> void allocate(T * &a, uint64_t n);

/*
 * Steps:
 * 1. Load the index, forward and backward
 * 2. Initialize the FPGA
 *
 */

int main(int argc, char *argv[]) {

    struct timeval  tv1, tv2;

    string ext = "ral";

    // variables for file IO
    //FILE * N_fp = NULL;
    FILE * FM_fp = NULL;
    FILE * FM_meta_fp = NULL;

    //FILE * N_fp_rev = NULL;
    FILE * FM_fp_rev = NULL;
    FILE * FM_meta_fp_rev = NULL;

    FILE *in_fp = NULL;
    FILE *out_fp = NULL;
    FILE *aligned_fp = NULL;

    uint64_t f_size;
	//uint64_t f_size_rev;

    // variables from meta-file
    uint64_t fmt_len;
    bool is32bit;

    uint32_t bucket_bwt_len;
    uint32_t bucket_bwt_len_rev;

    uint32_t bucket_pad_size;
    uint32_t bucket_pad_size_rev;

    uint64_t end_char_pos;
    uint64_t end_char_pos_rev;

    // variables for bwt/ FM-index
    uint32_t n_buckets;
    uint32_t n_buckets_rev;

    index32_t * index32;
    index32_t * index32_rev;

    uint32_t cnt32[FM_BP_RANGE + 1] 	= {0};
    uint32_t cnt32_rev[FM_BP_RANGE + 1] = {0};

    // number of chromosomes
    uint16_t chrs_num;
    uint16_t chrs_num_rev;
    
    // N char info
    uint64_t N_cluster;
    uint64_t N_cluster_rev;

#if IS_SIM == 0
    // program usage
    if (argc != 3) {
        printf("usage: %s <index basename> <Reads file>\n", argv[0]);
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

#elif IS_SIM == 1

    string argv_1 = "sample_index/sample922";
    string s1 = argv_1 + ".1." + ext;
    string s2 = argv_1 + ".2." + ext;
    string s3 = argv_1 + ".3." + ext;

    string s4 = argv_1 + ".4." + ext;
    string s5 = argv_1 + ".5." + ext;
    string s6 = argv_1 + ".6." + ext;

    string s7 = argv_1 + ".7." + ext;
    string s8 = argv_1 + ".8." + ext;

    string s9 = argv_1 + ".9." + ext;
    string s10 = argv_1 + ".10." + ext;

    

#endif

    string r1 = "om1w.fq";
    string ra1 = "om1w.txt";


    //** load index
    printf("Reading meta data ... \n"); fflush(stdout);

    openFile(&FM_fp, s1, "r");
    openFile(&FM_meta_fp, s2, "r");

    openFile(&FM_fp_rev, s4, "r");
    openFile(&FM_meta_fp_rev, s5, "r");

    f_size = fileSizeBytes(FM_fp);
    read_meta(FM_meta_fp, 
                &fmt_len, 
                &is32bit, 
                &bucket_bwt_len, 
                &end_char_pos, 
                &bucket_pad_size,
                &N_cluster,
                &chrs_num);

    //f_size_rev = fileSizeBytes(FM_fp_rev);
    read_meta(FM_meta_fp_rev, 
                &fmt_len, 
                &is32bit, 
                &bucket_bwt_len_rev, 
                &end_char_pos_rev, 
                &bucket_pad_size_rev,
                &N_cluster_rev,
                &chrs_num_rev);

    n_buckets 		= CEIL(fmt_len, bucket_bwt_len);
    n_buckets_rev 	= CEIL(fmt_len, bucket_bwt_len_rev);

    fclose(FM_meta_fp);
    fclose(FM_meta_fp_rev);

    printf("FINISH ---> Reading meta data\n\n");fflush(stdout);


    //** read index
    printf("Reading index ... "); fflush(stdout);
    gettimeofday(&tv1, NULL);

    if (is32bit){
        uint64_t tmp_cnt;
        //read the i-table (note it is always stored at 64bit)
        for (int i = 0; i < FM_BP_RANGE + 1; i++){
            readFile(FM_fp, & tmp_cnt, sizeof(uint64_t));
            cnt32[i] = (uint32_t) tmp_cnt; //cast it to 32bits

            readFile(FM_fp_rev, & tmp_cnt, sizeof(uint64_t));
            cnt32_rev[i]= (uint32_t) tmp_cnt;

        }
    }

    allocate(index32, n_buckets * BUCKET_SIZE);
    readFile(FM_fp, index32, n_buckets * BUCKET_SIZE);

    allocate(index32_rev, n_buckets_rev * BUCKET_SIZE);
    readFile(FM_fp_rev, index32_rev, n_buckets_rev * BUCKET_SIZE);

    fclose(FM_fp);
    fclose(FM_fp_rev);
    printf("FINISH ---> Reading index\n\n");fflush(stdout);


    //** program FPGA
    printf("\nInitializing FPGA ... \n"); fflush(stdout);

    // load one mismatch match kernel bitstream
    fpga_t fpga_var;
	fpga_var.maxfile = OM1Way_init();
	fpga_var.group = max_load_group(fpga_var.maxfile, MAXOS_EXCLUSIVE, "*", 1);
	fpga_var.engine = max_lock_any(fpga_var.group);

	printf("OK\n"); fflush(stdout);


    //** write index into FPGA
    printf("\nWriting both forward and backward index into FPGA... \n %d",n_buckets); fflush(stdout);
    std::thread writeindex_thread;
    writeindex_thread = std::thread(writeIndexBothWay, index32, index32_rev, n_buckets * BUCKET_SIZE, fpga_var);


	//writeIndex(index32, n_buckets * BUCKET_SIZE, 0, fpga_var);
	//writeIndex(index32_rev, n_buckets_rev * BUCKET_SIZE, n_buckets * BUCKET_SIZE, fpga_var);//#TODO: high_init_rev may not be 32bit

	printf("OK\n"); fflush(stdout);


	//** Loading the first batch of short reads
    printf("Loading short reads (1st batch) ... ");fflush(stdout);
    printf("Loading short reads (1st batch) ... \n");fflush(stdout);

    // allocate I/O buffers
    char * in_buff 	= new char [BUFF_SIZE + 512];
    char * out_buff = new char [BUFF_SIZE];
    char * aligned_buff = new char [BUFF_SIZE];

    std::vector<read2Bit_t> reads1, reads2;

    openFile(&out_fp, r1, "w+");
    openFile(&aligned_fp, ra1, "w+");

    // read first batch
#if IS_SIM == 0
    openFile(&in_fp, argv[2], "r");
#elif IS_SIM == 1
    openFile(&in_fp, "sample_index/112reads.fq", "r");
#endif

    f_size = fileSizeBytes(in_fp);

    // determine how much data to read from the file
    uint64_t bytes_r = 0;
    uint64_t size_r =  bytes_r + BUFF_SIZE <= f_size ? BUFF_SIZE : f_size - bytes_r;

    loadReads(in_fp, reads1, in_buff, size_r, &bytes_r, true);

    // print the first 5 reads
    for (uint32_t i = 0; i < 5; i++) {
        printf("%u: %s, %s\n", i, reads1[i].at_line.c_str(), reads1[i].seq);
    }
    printf("\n");

    // hopefully the write index thread will be finished by now
    printf("Waiting for index writing to finish ... \n"); fflush(stdout);
    writeindex_thread.join();
    printf("FINISH ---> Writing index\n\n"); fflush(stdout);

    printf("Aligning reads and Loading more short reads ... \n"); fflush(stdout);

	bmk_t bmk1, bmk2;
	bmk1.in_size = 0;
	bmk1.out_size = 0;
	bmk1.process_time = 0;
	bmk1.aligned_cnt = 0;

	bmk2.in_size = 0;
	bmk2.out_size = 0;
	bmk2.process_time = 0;
	bmk2.aligned_cnt = 0;


    uint64_t aligned_cnt1 = 0;
    uint64_t aligned_cnt2 = 0;
    uint32_t cnt = 0;

    std::thread loadread_thread;
    std::thread writeAligned_thread;

    for (char thread_cnt = 0; ; thread_cnt++) {

    	bool r_ctrl = bytes_r < f_size ? true : false;
    	size_r =  bytes_r + BUFF_SIZE <= f_size ? BUFF_SIZE : f_size - bytes_r;

    	// load to reads2, process reads1
        if (!(thread_cnt % 2)) {

        	printf("\nLoad Read (A) Again\n"); fflush(stdout);
        	loadread_thread = std::thread(loadReads,
        			in_fp,
					std::ref(reads2),
					in_buff,
					size_r,
					&bytes_r,
					true);

			 if (reads1.size() > 0) {
				cnt += reads1.size();

				gettimeofday(&tv1, NULL);

				align(reads1,
			            n_buckets * BUCKET_SIZE,//index_bytes
			            cnt32[FM_BP_RANGE],
						cnt32_rev[FM_BP_RANGE],
			            bucket_bwt_len,
						bucket_bwt_len_rev,
			            end_char_pos,
						end_char_pos_rev,
			            fpga_var, &bmk1);

                // writing aligned read to a plain text file for later processing
                writeAligned_thread = std::thread(writeAligned,
        			aligned_fp,
					std::ref(reads1),
					aligned_buff);

                aligned_cnt1 = aligned_cnt1 + writeReads(out_fp, reads1, out_buff);

                writeAligned_thread.join();

                gettimeofday(&tv2, NULL);
                printf("FINISH: One round alignmemt [%.2f s]\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
                        (double) (tv2.tv_sec - tv1.tv_sec));

                printf("...................One pass finished.......................................\n"); fflush(stdout);


			 }else{
				 break;
			 }
	    // load to reads1, process reads2
		}else{

			printf("\nLoad Read (B) Again\n"); fflush(stdout);
        	loadread_thread = std::thread(loadReads,
        			in_fp,
					std::ref(reads1),
					in_buff,
					size_r,
					&bytes_r,
					true);

			 if (reads2.size() > 0) {
				cnt += reads2.size();

				gettimeofday(&tv1, NULL);

				align(reads2,
			            n_buckets * BUCKET_SIZE,//index_bytes
			            cnt32[FM_BP_RANGE],
						cnt32_rev[FM_BP_RANGE],
			            bucket_bwt_len,
						bucket_bwt_len_rev,
			            end_char_pos,
						end_char_pos_rev,
			            fpga_var, &bmk2);

                // writing aligned read to a plain text file for later processing
                writeAligned_thread = std::thread(writeAligned,
        			aligned_fp,
					std::ref(reads2),
					aligned_buff);

				aligned_cnt2 = aligned_cnt2 + writeReads(out_fp, reads2, out_buff);

                writeAligned_thread.join();

				gettimeofday(&tv2, NULL);
        		printf("FINISH: One round alignmemt [%.2f s]\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
        					(double) (tv2.tv_sec - tv1.tv_sec));

        		printf("...................One pass finished.......................................\n"); fflush(stdout);

			 }else{
				 break;
			 }
		}
        loadread_thread.join();
    }
    loadread_thread.join();

    printf("processed %u reads\n", cnt);
    printf("aligned %lu reads (file write)\n", aligned_cnt1 + aligned_cnt2);
    printf("aligned %lu reads (alignment)\n", bmk2.aligned_cnt  + bmk1.aligned_cnt);
    printf("FPGA time in seconds: %.2f s\n", bmk2.process_time  + bmk1.process_time );

	max_unlock(fpga_var.engine);
  	max_unload_group(fpga_var.group);
  	max_file_free(fpga_var.maxfile);

	fclose(in_fp);
	fclose(out_fp);
    fclose(aligned_fp);

    delete [] index32;
    delete [] aligned_buff;

	return 0;
}

template<typename T> void allocate(T * &a, uint64_t n)
{
  a = new T [n/sizeof(T)];
  if (!a) {
	fprintf(stderr, "error: unable to allocate memory!\n");
    exit(1);
  }
}
