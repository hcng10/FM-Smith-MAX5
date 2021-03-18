#include <stdio.h>
#include <stdint.h>
#include <string>
#include <cstring>
#include <iostream>
#include <bitset>
#include <thread>

#include "fm-sa-def.h"
#include "readFM.h"
#include "reads.h"
#include "align.h"
#include "file_op.h"
#include "writeSAM.h"



using namespace std;

template<typename T> void allocate(T * &a, uint64_t n);

uint8_t charToBin(char s){
    switch(s) {
		case 'A': return 0; break;
		case 'C': return 1; break;
		case 'G': return 2; break;
		case 'T': return 3; break;
		case 'a': return 0; break;
		case 'c': return 1; break;
		case 'g': return 2; break;
		case 't': return 3; break;
    }
    return 0;

}

int main(int argc, char *argv[]) {

    struct timeval  tv1, tv2;
    string ext = "ral";

    bool is_seeding = true;

    // variables for file
    FILE * N_fp = NULL;
    FILE * FM_fp = NULL;
    FILE * FM_meta_fp = NULL;

    //FILE * N_fp_rev = NULL;
    //FILE * FM_fp_rev = NULL;
    //FILE * FM_meta_fp_rev = NULL;

    FILE * SA_fp = NULL;
    //FILE * rev_SA_fp = NULL;

    FILE * chr_fp = NULL;

    FILE * in_fp = NULL;
    FILE * out_fp = NULL;
    FILE * outN_fp = NULL;
    FILE * sam_fp = NULL;

    uint64_t f_size;

    // variables from meta-file
    uint64_t fmt_len;
    bool is32bit;

    uint32_t bucket_bwt_len;
    uint32_t bucket_pad_size;

    uint64_t end_char_pos;

    // variables for bwt/ FM-index
    uint32_t n_buckets;

    index32_t * index32;
    uint32_t cnt32[FM_BP_RANGE + 1] = {0};

    // suffix array
    uint32_t * sai;

    // number of chromosomes
    vector<chr_t> chrs;
    uint16_t chrs_num;
    
    // N char info
    //vector<nchar_cluster_t> nchar_clusters;
    uint64_t N_cluster;

#if IS_SIM == 0
    // program usage
    if (argc != 4) {
        printf("usage: %s <index basename> <Reads file> <is Seeding>\n", argv[0]);
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

    is_seeding = ((strcmp(argv[3], "true") == 0) | (string(argv[3]).compare("1") == 0)) ? true : false;

#elif IS_SIM == 1

    string argv_1 = "sample_index/sample_922_wN";
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
    string r1 = "em2w.fq";
    string r1N = "em2wN.fq";
    string rs1 = "em2w.txt";

    if (is_seeding == false){
        printf("** Alignment using FPGA with EM **\n");
    }
    else{
        printf("** Seeding using FPGA **\n");
    }

    fflush(stdout);

    // load index
    printf("Reading meta data ... "); fflush(stdout);

    openFile(&FM_fp, s1, "r");
    openFile(&FM_meta_fp, s2, "r");


    f_size = fileSizeBytes(FM_fp);
    read_meta(FM_meta_fp, 
                &fmt_len, 
                &is32bit, 
                &bucket_bwt_len,
                &end_char_pos, 
                &bucket_pad_size,
                &N_cluster,
                &chrs_num);


    n_buckets = CEIL(fmt_len, bucket_bwt_len);

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
    /*if (fread(sai, sizeof(uint32_t), fmt_len, SA_fp) != fmt_len) {
        fprintf(stderr, "error: unable to read SA file!\n");
        exit(1);
    }*/

    fclose(SA_fp);
    printf("FINISH ---> Reading SA\n\n");



    // read N char info
    //printf("Reading N char info in reference ... \n");fflush(stdout);

    //openFile(&N_fp, s3, "r");
    //readNinfo(N_fp, N_cluster, nchar_clusters);

    //fclose(N_fp);
    //printf("FINISH ---> Reading N char info\n\n");


    // read the name of chromosome
    printf("Reading the name of chromosome ... \n");fflush(stdout);

    openFile(&chr_fp, s9, "r");
    readChrName(chr_fp, chrs_num, chrs);

    fclose(chr_fp);
    printf("FINISH ---> Reading name of chromosome\n\n");



    printf("Reading index ... \n"); fflush(stdout);

    gettimeofday(&tv1, NULL);
    if (is32bit){
    	uint64_t tmp_cnt;
        //read the i-table (note it is always stored at 64bit)
        for (int i = 0; i < FM_BP_RANGE + 1; i++){

            readFile(FM_fp, & tmp_cnt, sizeof(uint64_t));
            cnt32[i] = (uint32_t) tmp_cnt; //cast it to 32bits
        }

        printf("\t Reading Buckets ... \n"); fflush(stdout);
        allocate(index32, n_buckets * BUCKET_SIZE);
        readFile(FM_fp, index32, n_buckets * BUCKET_SIZE);
    }
    gettimeofday(&tv2, NULL);
    printf("\t Reading index Time: [%.2f s]\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
            				(double) (tv2.tv_sec - tv1.tv_sec));

    fclose(FM_fp);
    printf("FINISH ---> Reading index\n\n");fflush(stdout);


    printf("Initializing FPGA ... \n"); fflush(stdout);
    

    // load exact match kernel bitstream
    fpga_t fpga_var;
    gettimeofday(&tv1, NULL);

	fpga_var.maxfile = EM2WaySeed_init();
	fpga_var.group = max_load_group(fpga_var.maxfile, MAXOS_EXCLUSIVE, "*", 1);
	fpga_var.engine = max_lock_any(fpga_var.group);

	gettimeofday(&tv2, NULL);

	printf("FINISH ---> Initializing FPGA [%.2f s]\n\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
			(double) (tv2.tv_sec - tv1.tv_sec)); fflush(stdout);


    //** write index into FPGA
    printf("Writing backward index into FPGA... \n"); fflush(stdout);
	max_run_t *exec_status = writeIndex(index32, n_buckets * BUCKET_SIZE, fpga_var);

	//std::thread writeIndex_thread;
	/*writeIndex_thread = std::thread(writeIndexBothWay,
									index32, index32_rev,
									n_buckets * BUCKET_SIZE,
									n_buckets,
									fpga_var);*/



    printf("\nLoading short reads (1st batch) ... \n");fflush(stdout);

    // allocate I/O buffers
    char * in_buff = new char [BUFF_SIZE + 512];
    char * out_buff = new char [BUFF_SIZE];
    char * sam_buff = new char [BUFF_SIZE];

    std::vector<read2Bit_t> reads1, reads2;
    std::vector<read3Bit_t> reads3, reads4;

    openFile(&out_fp, r1, "w+");
    openFile(&outN_fp, r1N, "w+");
    openFile(&sam_fp, rs1, "w+");

    // read first batch
#if IS_SIM == 0
    openFile(&in_fp, argv[2], "r");
#elif IS_SIM == 1
    openFile(&in_fp, "sample_index/reads_test_1.fq", "r");
#endif

    f_size = fileSizeBytes(in_fp);

    // determine how much data to read from the file
    uint64_t bytes_r = 0;
    uint64_t size_r =  bytes_r + BUFF_SIZE <= f_size ? BUFF_SIZE : f_size - bytes_r;

    loadReads(in_fp, reads1, reads3, in_buff, size_r, &bytes_r, true);

    for (uint32_t i = 0; i < 5; i++) {
        printf("%u: %s, %s\n", i, reads1[i].at_line.c_str(), reads1[i].seq);
    }
    printf("\n");

    if (reads3.size()) printf("Reads with 'N' characters\n");
    for (uint32_t i = 0; i < (reads3.size() < 5? reads3.size(): 5); i++) {
        printf("%u: %s, %s\n", i, reads3[i].at_line.c_str(), reads3[i].seq);
    }

    printf("Waiting for index writing to finish ... \n"); fflush(stdout);

    // check if the index is written into FPGA
    max_wait(exec_status);
    //writeIndex_thread.join();
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

    uint64_t N_cnt1 = 0;
    uint64_t N_cnt2 = 0;

    uint32_t cnt = 0;

    std::thread loadread_thread;
    std::thread writeSAM_thread;
    for (char thread_cnt = 0; ; thread_cnt++) {

    	bool r_ctrl = bytes_r < f_size ? true : false;
    	size_r =  bytes_r + BUFF_SIZE <= f_size ? BUFF_SIZE : f_size - bytes_r;

        // read to reads2, process reads1
        if (!(thread_cnt % 2)) {

        	printf("\nLoad Read (A) Again\n"); fflush(stdout);
        	loadread_thread = std::thread(loadReads,
        			in_fp,
					std::ref(reads2),
					std::ref(reads4),
					in_buff,
					size_r, &bytes_r, r_ctrl);

        	 if (reads1.size() > 0) {
                cnt += reads1.size();

                gettimeofday(&tv1, NULL);

                align(reads1,
                        n_buckets * BUCKET_SIZE,//index_bytes
                        cnt32[FM_BP_RANGE],
                        bucket_bwt_len,
                        end_char_pos,
                        is_seeding,
                        fpga_var,
                        &bmk1);

                // writing aligned read to SAM file
                /*writeSAM_thread = std::thread(convertSAM,
        			sam_fp,
					std::ref(reads1),
					sai,
					sam_buff,
					std::ref(chrs));*/
                if (is_seeding){
                    writeSAM_thread = std::thread(writeSeedPosDirectly_wThread,//writeSeedPosDirectly_wThread, writeSeedPosDirectly
                        sam_fp,
                        std::ref(reads1),
                        sai,
                        std::ref(chrs),
                        sam_buff
                    );
                }else{
                    writeSAM_thread = std::thread(writePosDirectly,//,
                        sam_fp,
                        std::ref(reads1),
                        sai,
                        std::ref(chrs),
                        sam_buff
                    );
                }


                aligned_cnt1 = aligned_cnt1 + writeReads(out_fp, reads1, is_seeding, out_buff);
                N_cnt1 = N_cnt1 + writeReadsN(outN_fp, reads3, out_buff);

                writeSAM_thread.join();

                gettimeofday(&tv2, NULL);
                printf("FINISH: One round alignmemt [%.2f s]\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
                        (double) (tv2.tv_sec - tv1.tv_sec));

                printf("...................One pass (A) finished.......................................\n"); fflush(stdout);


        	 }else
        		 break;

        }else{
 
        	printf("\nLoad Read (B) Again\n"); fflush(stdout);
        	loadread_thread = std::thread(loadReads,
        			in_fp,
					std::ref(reads1),
					std::ref(reads3),
					in_buff,
					size_r, &bytes_r, r_ctrl);


        	if (reads2.size() > 0) {
                cnt += reads2.size();

                gettimeofday(&tv1, NULL);

        		align(reads2,
        				n_buckets * BUCKET_SIZE,//index_bytes
						cnt32[FM_BP_RANGE],
						bucket_bwt_len,
						end_char_pos,
                        is_seeding,
						fpga_var,
						&bmk2);
                        
                // writing aligned read to SAM file
                /*writeSAM_thread = std::thread(convertSAM,
        			sam_fp,
					std::ref(reads2),
					sai,
					sam_buff,
					std::ref(chrs));*/

                if (is_seeding){
                    writeSAM_thread = std::thread(writeSeedPosDirectly_wThread,//writeSeedPosDirectly_wThread,
                        sam_fp,
                        std::ref(reads2),
                        sai,
                        std::ref(chrs),
                        sam_buff
                    );
                }else{
                     writeSAM_thread = std::thread(writePosDirectly,//,
                        sam_fp,
                        std::ref(reads2),
                        sai,
                        std::ref(chrs),
                        sam_buff
                    );                   
                }

        		aligned_cnt2 = aligned_cnt2 + writeReads(out_fp, reads2, is_seeding, out_buff);
        		N_cnt2 = N_cnt2 + writeReadsN(outN_fp, reads4, out_buff);
                
                writeSAM_thread.join();

                gettimeofday(&tv2, NULL);
                printf("FINISH: One round alignmemt [%.2f s]\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
                            (double) (tv2.tv_sec - tv1.tv_sec));

                printf("...................One (B) pass finished.......................................\n"); fflush(stdout);

        	}else
        		break;

        }
        loadread_thread.join();


    }
    loadread_thread.join();


    printf("processed %lu reads\n", cnt + N_cnt1  + N_cnt2);
    printf("aligned %lu reads (file write)\n", aligned_cnt1 + aligned_cnt2);
    printf("aligned %lu reads (alignment)\n", bmk2.aligned_cnt  + bmk1.aligned_cnt);


    printf("unaligned %lu N char reads \n", N_cnt1  + N_cnt2);
    printf("FPGA time in seconds: %.2f s\n", bmk2.process_time  + bmk1.process_time );

    fflush(stdout);


	//cleanup
	max_unlock(fpga_var.engine);
	max_unload_group(fpga_var.group);
	max_file_free(fpga_var.maxfile);

    printf("Unlock DFE \n");fflush(stdout);

	fclose(in_fp);
	fclose(out_fp);
	fclose(outN_fp);
    fclose(sam_fp);

    printf("File closed\n");fflush(stdout);
    
    //delete [] in_buff;cerr<<"chk1\n";
    //delete [] out_buff;
    delete [] sam_buff;
    //delete [] index32;

    delete [] sai;


    //align(reads1, cnt32, index32, n_buckets * BUCKET_SIZE, cnt32[FM_BP_RANGE], bucket_bwt_len, end_char_pos);


}

template<typename T> void allocate(T * &a, uint64_t n)
{
  a = new T [n/sizeof(T)];
  if (!a) {
	fprintf(stderr, "error: unable to allocate memory!\n");
    exit(1);
  }
}
