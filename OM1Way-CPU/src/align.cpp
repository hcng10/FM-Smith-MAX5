/*
 * align.cpp
 *
 *  Created on: Feb 6, 2019
 *      Author: hn915
 */


#include "align.h"


/*
 * Function declaration
 */
// get partition sizes
void getPartSize(uint32_t *part_size, uint32_t items);
//inline char getPckVal(uint8_t pck_val, uint8_t idx);
//inline uint8_t getPckQS(uint16_t pck_qs, uint8_t idx);

void actionWrite(index32_t *index, uint64_t index_bytes, uint64_t offset, max_engine_t *engine);
// kernel read action
void actionRead(out_t **out, uint64_t *hit_count, uint64_t offset,
		max_engine_t *engine);

void actionAlign(in_t **in,
					uint32_t *part_size,
					uint64_t index_bytes,
					uint32_t high_init,
					uint32_t high_init_rev,
					uint32_t end_char_bucket,
					uint32_t end_char_bucketi,
					uint32_t end_char_bucket_rev,
					uint32_t end_char_bucketi_rev,
					out_t **out,
					bool isReverseCmpt,
					max_engine_t *engine, bmk_t * bmk);

/*
 * Function definition
 */
/*inline char getPckVal(uint8_t pck_val, uint8_t idx){
	char sym;
	//uint8_t tmp = pck_val >> ((idx % (8/FM_BP_BIT))* FM_BP_BIT);
	uint8_t tmp = pck_val >> ((idx * FM_BP_BIT) % (sizeof(uint8_t)*8));
	uint8_t symVal = tmp & 3;

    switch(symVal) {
    case 0: sym = 'A'; break;
    case 1: sym = 'C'; break;
    case 2: sym = 'G'; break;
    case 3: sym = 'T'; break;
    default : sym = 'A';
    }
    return sym;
}


inline uint8_t getPckQS(uint16_t pck_qs, uint8_t idx){
	uint8_t qs;
	uint16_t tmp = pck_qs >> ((idx * FM_QS_BIT) % (sizeof(uint16_t)*8));
	qs = (uint8_t) (tmp & 63);

	return qs;
}*/

inline char getPckSym(uint8_t pck_val){
	char sym;
	uint8_t tmp = pck_val >> 6;
	tmp = tmp & 3;

    switch(tmp) {
    case 0: sym = 'A'; break;
    case 1: sym = 'C'; break;
    case 2: sym = 'G'; break;
    case 3: sym = 'T'; break;
    default : sym = 'A';
    }
    return sym;

}


// get partition size
void getPartSize(uint32_t *part_size, uint32_t item_size){

	for (uint8_t i = 0; i < N_KRNL; i++) {
		if (item_size % (N_KRNL)!=0) {
			part_size[i] = CEIL(item_size, N_KRNL);
			item_size -= 1;
		}
		else
			part_size[i] = item_size/N_KRNL;
	}
}

void writeIndexBothWay(index32_t * index, index32_t * index_rev, uint64_t index_bytes, fpga_t fpga_var){
	struct timeval  tv1, tv2;

	// write index to DRAM memory
	printf("\t Transfer index (both ways) to DRAM ... \n"); fflush(stdout);
	gettimeofday(&tv1, NULL);
	/*****FPGA *************************************************/
	actionWrite(index, index_bytes, 0, fpga_var.engine);
	actionWrite(index_rev, index_bytes, index_bytes, fpga_var.engine);

	gettimeofday(&tv2, NULL);
	printf("FINISH ---> Transfer index [%.2f s]\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
			(double) (tv2.tv_sec - tv1.tv_sec)); fflush(stdout);

	return;

}

void writeIndex(index32_t * index, uint64_t index_bytes, uint64_t offset, fpga_t fpga_var){

	struct timeval  tv1, tv2;

	// write index to DRAM memory
	printf("\ttransferring index to DRAM ... "); fflush(stdout);
	gettimeofday(&tv1, NULL);
	/*****FPGA *************************************************/
	actionWrite(index, index_bytes, offset, fpga_var.engine);
	//TODO: add a reverse index

	gettimeofday(&tv2, NULL);
	printf("OK [%.2f s]\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
			(double) (tv2.tv_sec - tv1.tv_sec));

	return;
}

void parse_result(uint32_t * part_size,
				std::vector<read2Bit_t> &reads,
                out_t **out,
				bmk_t * bmk){

	struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);
    // get the number of threads to allow parallel processing when data is moved from
	// reads to part_size
	int n_threads = omp_get_max_threads();

    for (uint8_t i = 0; i < N_KRNL; i++){
#pragma omp parallel for num_threads(n_threads)
	    for (uint32_t j = 0; j < part_size[i]; j++) {

	        hit_t tmp;

	        uint32_t id = out[i][j].id;
	    	uint8_t n_hits = out[i][j].n_hits;

			if (n_hits == 0) continue;

			tmp.n_hits = n_hits;
	    	memcpy(tmp.mis_pos_sorted, out[i][j].mis_pos_sorted, sizeof(tmp.mis_pos_sorted));
	    	memcpy(tmp.low_sorted, out[i][j].low_sorted, sizeof(tmp.low_sorted));
	    	memcpy(tmp.high_sorted, out[i][j].high_sorted, sizeof(tmp.high_sorted));

	    	for (uint8_t k = 0; k < n_hits; k++){
	    	    uint8_t tmp_qs_mis_sym = out[i][j].qs_mis_sym_sorted[N_HITS-1-k];
	    	    tmp.mis_sym_sorted[N_HITS-1-k] = getPckSym(tmp_qs_mis_sym);
	    	    tmp.qs_sorted[N_HITS-1-k] = (tmp_qs_mis_sym & 63) +33;

	    	    //printf("Num: %d, ID: %d, Missym %c, QS %d \n", k, id, tmp.mis_sym_sorted[N_HITS-1-k], tmp.qs_sorted[N_HITS-1-k]);
	    	    //printf("Missym %c, QS %d \n", tmp.mis_sym_sorted[N_HITS-1-k], tmp.qs_sorted[N_HITS-1-k]);

	    	}
			//cout<<j<<" id: "<<id<<" n_hits "<<(int)n_hits<<" "<<(int)tmp.mis_pos_sorted[0]<<
				//" "<<(int)tmp.mis_pos_sorted[31]<<"\n"<<
				//(int)tmp.low_sorted[31]<<" "<<(int)tmp.high_sorted[31]<<" "<<
				//(char)tmp.mis_sym_sorted[31]<<" "<<(char)(tmp.qs_sorted[31])<<"\n";

	    	reads[id].n_hits = reads[id].n_hits + n_hits;
	    	reads[id].is_align = n_hits > 0? true: false;

	    	if (n_hits > 0){
	    		reads[id].hits.push_back(tmp);
	    	}

	    }
    }

	uint32_t total_aligned = 0;
	uint32_t total_unaligned = 0;

	uint32_t total_cnt = 0;

	for (uint8_t i = 0; i < N_KRNL; i++) {
		for (uint32_t j = 0; j < part_size[i]; j++) {

			total_cnt = total_cnt + 1;

			if (out[i][j].n_hits > 0){
				total_aligned = total_aligned + 1;
			}
			else{
				total_unaligned = total_unaligned + 1;
			}
		}
	}

	bmk->aligned_cnt = bmk->aligned_cnt + total_aligned;

	gettimeofday(&tv2, NULL);
	printf("OK [%.2f s], Total Aligned: %d, Total UnAligned %d total %d \n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
	 (double) (tv2.tv_sec - tv1.tv_sec),
	 total_aligned,
	 total_unaligned,
	 total_cnt);
}

void align(std::vector<read2Bit_t> &reads,
			//uint32_t * cnt32,
			//uint32_t * cnt32_rev,
			uint64_t index_bytes,
			//uint64_t index_bytes_rev,
			uint32_t high_init,
			uint32_t high_init_rev,
			uint32_t bucket_bwt_len,
			uint32_t bucket_bwt_len_rev,
			uint32_t end_char_pos,
			uint32_t end_char_pos_rev,
			fpga_t fpga_var,
			bmk_t * bmk){


	//struct timeval  tv1, tv2;

	in_t *in[N_KRNL];
	out_t *out1[N_KRNL];
	out_t *out2[N_KRNL];

	uint32_t part_size[N_KRNL];
	uint64_t hit_count[N_KRNL];

	uint32_t end_char_bucket	= (uint32_t) (end_char_pos / bucket_bwt_len);
	uint32_t end_char_bucketi  	= (uint32_t) (end_char_pos % bucket_bwt_len);

	uint32_t end_char_bucket_rev 	= (uint32_t) (end_char_pos_rev / bucket_bwt_len_rev);
	uint32_t end_char_bucketi_rev = (uint32_t) (end_char_pos_rev % bucket_bwt_len_rev);


    // get the number of threads to allow parallel processing when data is moved from
	// reads to part_size
	int n_threads = omp_get_max_threads();
	std::thread parse_thread;

	// generate kernel input by dividing them into smaller copies based on the number of kernels
	// size refers to the number of elements in the vector
	printf("\t splitting data into parts ... "); fflush(stdout);
	//gettimeofday(&tv1, NULL);
	getPartSize(part_size, reads.size());


	uint32_t rds_offset = 0;
	for (uint8_t i = 0; i < N_KRNL; i++) {
		//in[i] = new in_t [part_size[i]]; //+ latency];
		in[i] = new in_t [part_size[i] + LATENCY];

		if (!in[i]) {
			fprintf(stderr, "error: unable to allocate memory!\n");
			exit(1);
		}
		memset(in[i], 0, part_size[i] + LATENCY * sizeof(in_t));

#pragma omp parallel for num_threads(n_threads)
		for (uint32_t j = 0; j < part_size[i]; j++){
			uint32_t id = j + rds_offset;
			in[i][j].id = id;
			in[i][j].len = reads[id].seq_len;
			memcpy(in[i][j].pck_sym_qs, reads[id].pck_sym_qs, reads[id].seq_len);
			//memcpy(in[i][j].pck_sym, reads[id].pck_sym, CEIL(reads[id].seq_len, FM_BP_RANGE));
			in[i][j].is_pad = 0;

			//cout<<"in[i][j].id "<<j<<" "<<in[i][j].id<<" "<<(int)in[i][j].len<<"\n";
		}

		for (uint32_t j = 0; j < LATENCY; j++) {
			if (j < LATENCY - 1){
				// with 1, that means you wouldn't issue any memory request,
				// and hence when you see is_pad == 2, there should be no more
				// memory request
				in[i][part_size[i]+j].id = rds_offset + part_size[i] + j;
				in[i][part_size[i]+j].is_pad = 1;
			}else{
				in[i][part_size[i]+j].id = rds_offset + part_size[i] + j;
				in[i][part_size[i]+j].is_pad = 2;
			}

		}
		rds_offset += part_size[i];
	}
	printf("OK");

	// one mismatch align reads
	printf("\t aligning reads for 1 mismatch ... "); fflush(stdout);
	//gettimeofday(&tv1, NULL);
	memset(hit_count, 0, N_KRNL* sizeof(uint64_t));


	// align reads (FPGA)
	//gettimeofday(&tv1, NULL);
	for (uint8_t i = 0; i < N_KRNL; i++){
		out1[i] = new out_t [part_size[i]];
		out2[i] = new out_t [part_size[i]];
	    if (!out1[i] || !out2[i]) {
	    	 fprintf(stderr, "error: unable to allocate memory!\n");
	    	 exit(1);
	    }
	    memset(out1[i], 0, part_size[i] * sizeof(out_t));
	    memset(out2[i], 0, part_size[i] * sizeof(out_t));
	}

	actionAlign(in,
				part_size,
				index_bytes,
				high_init,
				high_init_rev,
				end_char_bucket,
				end_char_bucketi,
				end_char_bucket_rev,
				end_char_bucketi_rev,
				out1, false, fpga_var.engine, bmk);

	printf("OK"); fflush(stdout);

	// parse output
	printf("\tparsing results and alignment from forward direction... "); fflush(stdout);
parse_result(part_size, reads, out1, bmk);
	//parse_thread = std::thread(parse_result, part_size, std::ref(reads), out1, std::ref(bmk));

	actionAlign(in,
				part_size,
				index_bytes,
				high_init,
				high_init_rev,
				end_char_bucket,
				end_char_bucketi,
				end_char_bucket_rev,
				end_char_bucketi_rev,
				out2, true, fpga_var.engine, bmk);

	//parse_thread.join();

	printf("\tparsing final results ... \n"); fflush(stdout);

	parse_result(part_size, reads, out2, bmk);


	// cleanup
	for (int i = 0; i < N_KRNL; i++) {
		delete[] in[i];
		delete[] out1[i];
		delete[] out2[i];
	}
}


void actionWrite(index32_t *index, uint64_t index_bytes, uint64_t offset, max_engine_t *engine){

	OM1Way_Write_actions_t *om_write = new OM1Way_Write_actions_t;

	om_write->param_nBytes = index_bytes;
	om_write->param_offset = offset;
	// index is based on index32_t, so its index has to be divided by BURST_BYTES
	//om_write->instream_OMIndexToMger = (uint64_t*) &index[offset / BURST_BYTES];
	om_write->instream_OMIndexToMger = (uint64_t*) index;

	OM1Way_Write_run(engine, om_write);
	delete om_write;

}
//TODO trim down the unused variables
void actionAlign(in_t **in,
					uint32_t *part_size,
					uint64_t index_bytes,
					uint32_t high_init,
					uint32_t high_init_rev,
					uint32_t end_char_bucket,
					uint32_t end_char_bucketi,
					uint32_t end_char_bucket_rev,
					uint32_t end_char_bucketi_rev,
					out_t **out,
					bool isReverseCmpt,
					max_engine_t *engine, bmk_t * bmk){

	 struct timeval  tv1, tv2;
	 double time_taken;


	  uint64_t n_bytes_in[N_KRNL];
	  uint64_t n_bytes_out[N_KRNL];

	  for (uint8_t i = 0; i < N_KRNL; i++) {
		  n_bytes_in[i] = (part_size[i] + LATENCY) * sizeof(in_t);
	      n_bytes_out[i] = part_size[i] * sizeof(out_t);//TODO: supposed to be changed
	  }



	  OM1Way_Align_actions_t * om_align = new OM1Way_Align_actions_t;

	  om_align->param_endCharBucket = end_char_bucket;
	  om_align->param_endCharBucketi = end_char_bucketi;

	  om_align->param_endCharBucketRev = end_char_bucket_rev;
	  om_align->param_endCharBucketiRev = end_char_bucketi_rev;

	  om_align->param_highInit = high_init;
	  om_align->param_highInitRev = high_init_rev;

	  om_align->param_offset = LATENCY;

	  //om_align->param_phase = phase;

	  //index is based on memory burst, so its index has to be divided by BURST_BYTES
	  om_align->param_readOffset = index_bytes / BURST_BYTES;
	  om_align->param_isReverseCmpt = (uint8_t)isReverseCmpt;
	  	  	  	  	  /*om_align->param_writeOffset = (2 * index_bytes) / BURST_BYTES+10;*/

	  om_align->param_nBytesInput = n_bytes_in;
	  om_align->param_nBytesOutput = n_bytes_out;

#if N_KRNL > 0
	  om_align->instream_readIn0 = (uint64_t*)in[0];
	  	  	  	  	  /*om_align->outscalar_OMKernel0_hitCount = & hit_count[0];*/
	  om_align->outstream_omOut0 = (uint64_t*)out[0];
#endif

#if N_KRNL > 1
	  om_align->instream_readIn1 = (uint64_t*)in[1];
	  om_align->outstream_omOut1 = (uint64_t*)out[1];
#endif

#if N_KRNL > 2
	  om_align->instream_readIn2 = (uint64_t*)in[2];
	  om_align->outstream_omOut2 = (uint64_t*)out[2];
#endif

	gettimeofday(&tv1, NULL);
	OM1Way_Align_run(engine, om_align);
    gettimeofday(&tv2, NULL);

    time_taken = ((double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
          (double) (tv2.tv_sec - tv1.tv_sec));

    printf("**FPGA TIME:** [%.2f s]\n", time_taken);
    fflush(stdout);

    bmk->process_time = bmk->process_time + time_taken;

	delete om_align;

}



// kernel read action
/*void actionRead(out_t **out, uint64_t *hit_count, uint64_t offset,
		max_engine_t *engine){

	uint64_t n_bytes[N_KRNL];
	for (uint8_t i = 0; i < N_KRNL; i++) {
		n_bytes[i] = CEIL(hit_count[i], N_HITS) * N_HITS * sizeof(out_t);
	}

	OMmk1_Read_actions_t *om_read = new OMmk1_Read_actions_t;
	om_read->param_nBytes = n_bytes;
	om_read->param_offset = offset;


#if N_KRNL > 0
    om_read->outstream_outputToHost0 = (uint64_t*)out[0];
#endif

#if N_KRNL > 1
    om_read->outstream_outputToHost1 = (uint64_t*)out[1];
#endif

#if N_KRNL > 2
    om_read->outstream_outputToHost1 = (uint64_t*)out[2];
#endif

    OMmk1_Read_run(engine, om_read);
    delete om_read;

}*/

