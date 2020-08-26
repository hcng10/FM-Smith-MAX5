/*
 * align.cpp
 *
 *  Created on: Feb 6, 2019
 *      Author: hn915
 */


#include "align.h"



// get partition sizes
void getPartSize(uint32_t *part_size, uint32_t items);

max_run_t * actionWrite(index32_t *index, uint64_t index_bytes, uint64_t offset, max_engine_t *engine);
void actionAlign(in_t **in,
					out_t **out,
					uint32_t * part_size,
					//uint64_t index_bytes,
					uint32_t high_init,
					//uint32_t high_init_rev,
					uint32_t end_char_bucket,
					uint32_t end_char_bucketi,
					//uint32_t end_char_bucket_rev,
					//uint32_t end_char_bucketi_rev,
					max_engine_t * engine, bmk_t * bmk);



/*void writeIndexBothWay(index32_t * index32, index32_t * index32_rev,
						uint64_t index_bytes,
						uint32_t n_buckets,
						fpga_t fpga_var){

	printf("THREAD: Writing backward index ... \n"); fflush(stdout);
	actionWrite(index32, index_bytes, 0, fpga_var.engine);
	//printf("THREAD: Writing forward index ... \n"); fflush(stdout);
	actionWrite(index32_rev, index_bytes, n_buckets * BUCKET_SIZE, fpga_var.engine);
}*/


max_run_t *writeIndex(index32_t * index, uint64_t index_bytes, fpga_t fpga_var){

	// write index to DRAM memory
	printf("\t Transfer index to DRAM ... \n"); fflush(stdout);
	/*****FPGA *************************************************/
	max_run_t *exec_status = actionWrite(index, index_bytes, 0, fpga_var.engine);

	return exec_status;
}


void align(std::vector<read2Bit_t> &reads,
		uint64_t index_bytes,
		uint32_t high_init,
		uint32_t bucket_bwt_len,
		uint32_t end_char_pos,
		fpga_t fpga_var,
		bmk_t * bmk){

	struct timeval  tv1, tv2;
	//gettimeofday(&tv1, NULL);

	in_t *in[N_KRNL];
	out_t *out[N_KRNL];

	uint32_t part_size[N_KRNL];

	uint32_t end_char_bucket = (uint32_t) (end_char_pos / bucket_bwt_len);
	uint32_t end_char_bucketi = (uint32_t) (end_char_pos % bucket_bwt_len);
	

	// get the number of threads to allow parallel processing when data is moved from
	// reads to part_size
	int n_threads = omp_get_max_threads();


	// generate kernel input by dividing them into smaller copies based on the number of kernels
	// size refers to the number of elements in the vector
	printf("\t splitting data into parts ..."); fflush(stdout);
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
		//memset(in[i], 0, part_size[i]);
		memset(in[i], 0, part_size[i] + LATENCY * sizeof(in_t));

#pragma omp parallel for num_threads(n_threads)
		for (uint32_t j = 0; j < part_size[i]; j++){
			uint32_t id = j + rds_offset;
			in[i][j].id = id;
			in[i][j].len = reads[id].seq_len;
			memcpy(in[i][j].pck_sym, reads[id].pck_sym, CEIL(reads[id].seq_len, FM_BP_RANGE));
			in[i][j].is_pad = 0;

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

	printf("OK\n");
	fflush(stdout);

	// exact align reads (FPGA)
	printf("\t exact aligning reads ... "); fflush(stdout);


	for (uint8_t i = 0; i < N_KRNL; i++){
		out[i] = new out_t [part_size[i]];
	    if (!out[i]) {
	    	 fprintf(stderr, "error: unable to allocate memory!\n");
	    	 exit(1);
	    }
	    memset(out[i], 0, part_size[i] * sizeof(out_t));
	}

	actionAlign(in,
				out,
				part_size,
				high_init,
				end_char_bucket,
				end_char_bucketi,
				fpga_var.engine, bmk);
				//index_bytes,
				//high_init_rev,
				//end_char_bucket_rev,
				//end_char_bucketi_rev,

	printf("OK\n");

	// parse output
	printf("\t parsing results ... "); fflush(stdout);
	gettimeofday(&tv1, NULL);

	for (uint8_t i = 0; i < N_KRNL; i++) {
#pragma omp parallel for num_threads(n_threads)
		for (uint32_t j = 0; j < part_size[i]; j++) {
			uint32_t id = out[i][j].id;

			reads[id].low_bw = out[i][j].low_bw;
			reads[id].high_bw = out[i][j].high_bw;

			reads[id].low_fw = out[i][j].low_fw;
			reads[id].high_fw = out[i][j].high_fw;


			if (out[i][j].isaligned_bw == 1 ){
				reads[id].isaligned_bw = true;
			}else{
				reads[id].isaligned_bw = false;
			}

			if (out[i][j].isaligned_fw == 1 ){
				reads[id].isaligned_fw = true;
			}else{
				reads[id].isaligned_fw = false;
			}

			//cout<<"RESULT: "<<id<<" low_bw: "<<out[i][j].low_bw<<" high_bw: "<<out[i][j].high_bw<<" low_fw: "<<out[i][j].low_fw<<" high_fw: "<<out[i][j].high_fw<<"\n";
			//cout<<"RESULT: "<<id<<" isaligned_bw " <<reads[id].isaligned_bw<<" isaligned_fw "<<reads[id].isaligned_fw<<"\n";
		}
	}

	uint32_t total_aligned = 0;
	uint32_t total_unaligned = 0;

	uint32_t total_cnt = 0;



	for (uint8_t i = 0; i < N_KRNL; i++) {
		for (uint32_t j = 0; j < part_size[i]; j++) {

			total_cnt = total_cnt + 1;

			if ((out[i][j].isaligned_bw == 1) || (out[i][j].isaligned_fw == 1)){
				total_aligned = total_aligned + 1;
			}
			else{
				total_unaligned = total_unaligned + 1;
			}
		}
	}


	bmk->aligned_cnt = bmk->aligned_cnt + total_aligned;

	gettimeofday(&tv2, NULL);
	printf("[%.2f s] OK, Total Aligned: %d, Total UnAligned %d total %d \n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
	 (double) (tv2.tv_sec - tv1.tv_sec),
	 total_aligned,
	 total_unaligned,
	 total_cnt);


	for (int i = 0; i < N_KRNL; i++) {
		delete[] in[i];
		delete[] out[i];
	}

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

max_run_t * actionWrite(index32_t *index, uint64_t index_bytes, uint64_t offset, max_engine_t *engine){

	EM2Way_Write_actions_t *em_write = new EM2Way_Write_actions_t;

	em_write->param_nBytes = index_bytes;
	em_write->param_offset = offset;
	em_write->instream_EMIndexToMger = (uint64_t*)index;


	max_run_t *exec_status = EM2Way_Write_run_nonblock(engine, em_write);
	//EM2Way_Write_run(engine, em_write);
	delete em_write;

	return exec_status;
}

void actionAlign(in_t **in,
					out_t **out,
					uint32_t * part_size,
					uint32_t high_init,
					uint32_t end_char_bucket,
					uint32_t end_char_bucketi,
					max_engine_t * engine, bmk_t * bmk){
					//uint64_t index_bytes,
					//uint32_t high_init_rev,
					//uint32_t end_char_bucket_rev,
					//uint32_t end_char_bucketi_rev,

	struct timeval  tv1, tv2;
	double time_taken;
	uint64_t bytes_in_size = 0;
	uint64_t bytes_out_size = 0;

	// determine how many bytes we need for input and output
	uint64_t n_bytes_in[N_KRNL];
	uint64_t n_bytes_out[N_KRNL];

	// get the number of bytes from the array parts
	for (uint8_t i = 0; i < N_KRNL; i++) {
		n_bytes_in[i] = (part_size[i] + LATENCY) * sizeof(in_t);
		n_bytes_out[i] = part_size[i] * sizeof(out_t);

		bytes_in_size = bytes_in_size + n_bytes_in[i];
		bytes_out_size = bytes_out_size + n_bytes_out[i];

	}

	EM2Way_Align_actions_t *em_align = new EM2Way_Align_actions_t;

	em_align->param_offset = LATENCY;

	em_align->param_highInit = high_init;

	em_align->param_endCharBucket = end_char_bucket;
	em_align->param_endCharBucketi = end_char_bucketi;


    em_align->param_nBytesInput = n_bytes_in;
    em_align->param_nBytesOutput = n_bytes_out;


    //printf("\n %ld the number of bytes %d %d\n", n_bytes_in[0], end_char_bucket, end_char_bucketi); fflush(stdout);

    em_align->instream_readIn0 = (uint64_t*)in[0];
    em_align->outstream_alignOut0 = (uint64_t*)out[0];

#if N_KRNL > 1
    em_align->instream_readIn1 = (uint64_t*)in[1];
    em_align->outstream_alignOut1 = (uint64_t*)out[1];
#endif

#if N_KRNL > 2
    em_align->instream_readIn2 = (uint64_t*)in[2];
    em_align->outstream_alignOut2 = (uint64_t*)out[2];
#endif

    gettimeofday(&tv1, NULL);
    EM2Way_Align_run(engine, em_align);
    gettimeofday(&tv2, NULL);

    time_taken = ((double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
          (double) (tv2.tv_sec - tv1.tv_sec));

    cout<<"[**FPGA TIME:** "<< time_taken<<"] "; //<< " byte in: "<< bytes_in_size << " byte out: " << bytes_out_size <<"\n";
    fflush(stdout);

    bmk->in_size =  bmk->in_size + bytes_in_size;
    bmk->out_size =  bmk->out_size + bytes_out_size;
    bmk->process_time = bmk->process_time + time_taken;

    delete em_align;

}

