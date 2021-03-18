
#ifndef ALIGN_H_
#define ALIGN_H_

#include <math.h>
#include "fm-sa-def.h"
#include "Maxfiles.h"



struct bmk_t{
	uint64_t in_size;
	uint64_t out_size;
	double process_time;
	uint64_t aligned_cnt;
};

struct fpga_t{
	max_file_t * maxfile;
	max_engine_t * engine;
	max_group_t * group;
};

// kernel input
struct in_t {
	uint32_t id;
	uint8_t read[CEIL(READ_DEF_LEN * BP_BIT, 8)];
	uint8_t readLen;
	uint8_t ref[CEIL(REF_DEF_LEN * BP_BIT, 8)];
	uint8_t refLen;
	uint8_t qs[CEIL(READ_DEF_LEN * QS_BIT, 8)];
	uint8_t isPad[PADSIZE/8];
};


// kernel output
struct out_t{
	uint32_t id;
	uint8_t rsltType [CEIL((READ_DEF_LEN + REF_DEF_LEN - 1 + RSLT_TAIL_ALIGN) * 2, 8)];
	uint16_t rsltMaxScore;
	uint32_t rsltLoc;
	uint8_t isOvrThld;
	uint8_t isPad[RSLT_PADSIZE/8];
};


void actionAlign(in_t **in,
					uint32_t *part_size,
					out_t **out,
					max_engine_t *engine, 
					bmk_t * bmk);


// get partition size
void getPartSize(uint32_t *part_size, uint64_t total_extend_cnt){

	uint64_t item_size = total_extend_cnt;

	for (uint8_t i = 0; i < N_KRNL; i++) {
		if (item_size % (N_KRNL)!=0) {
			part_size[i] = CEIL(item_size, N_KRNL);
			item_size -= 1;
		}
		else
			part_size[i] = item_size/N_KRNL;
	}
}

void parse_result(std::vector<Read_t> &reads, 
			std::vector<InToRead_t> & intoread, 
			uint32_t * part_size,
			out_t **out,
			uint32_t & total_aligned 
			){
	

	for (uint8_t i = 0; i < N_KRNL; i++){
	    for (uint32_t j = 0; j < part_size[i]; j++) {

			uint32_t id = out[i][j].id;
			bool is_aligned = (bool) out[i][j].isOvrThld;
			uint16_t aligned_score = out[i][j].rsltMaxScore;
			
			if (is_aligned == true){
				uint64_t read_idx = intoread[id].read_idx;
				uint16_t seed_idx = intoread[id].seed_idx;
				uint16_t extend_idx = intoread[id].extend_idx;

				Read_t & read = reads[read_idx];

				// record meta info
				if (read.q.empty() == true){
					total_aligned = total_aligned + 1;
				}

				AlignedInfo_t tmp_ainfo;
				tmp_ainfo.seed_idx = seed_idx;
				tmp_ainfo.extend_idx = extend_idx;
				tmp_ainfo.aligned_score = aligned_score;

				read.q.push(tmp_ainfo);
				
			}
//#if IS_SIM == 1
			//cout<<read_idx<<" all ID "<<id<<" "<<is_aligned<<" aligned score: "<<aligned_score<<" \n";
//#endif
		}
	}

	/*for (int i = 0; i < reads.size();i++){
		cout<<i<<" "<<reads[i].q.size();

		int tmp_size =  reads[i].q.size();

		for (uint k = 0;k<tmp_size;k++){
			AlignedInfo_t tmp = reads[i].q.top();
			
			cout<<" "<<(int)tmp.aligned_score;
			reads[i].q.pop();
		}
		cout<<"\n";
	}*/
}


void align(std::vector<Read_t> &reads, 
			std::vector<InToRead_t> & intoread, 
			uint64_t total_extend_cnt,
			fpga_t fpga_var,
			bmk_t * bmk){

	struct timeval  tv1, tv2;

	in_t *in[N_KRNL];
	out_t *out[N_KRNL];

	uint32_t part_size[N_KRNL];

	// get the number of threads to allow parallel processing when data is moved from
	// reads to in
	int n_threads = omp_get_max_threads();
	std::thread parse_thread;

	printf("\t splitting data into parts ... "); fflush(stdout);
	//gettimeofday(&tv1, NULL);
	getPartSize(part_size, total_extend_cnt);
//cout<<"************total extend out "<<total_extend_cnt<<"\n";

	uint32_t rds_offset = 0;
	for (uint8_t i = 0; i < N_KRNL; i++) {

		in[i] = new in_t [part_size[i] + LATENCY*2];//hack
		
		if (!in[i]) {
			fprintf(stderr, "error: unable to allocate memory!\n");
			exit(1);
		}
		memset(in[i], 0, (part_size[i] + LATENCY*2) * sizeof(in_t));


#pragma omp parallel for num_threads(n_threads)
		for (uint32_t j = 0; j < part_size[i]; j++){
			
			// id to get from vector InToRead_t
			uint32_t id = j + rds_offset;
			in[i][j].id = id;

			// based on InToRead_t, we can get back the read and the seed
			uint64_t read_idx;
			uint16_t seed_idx;
			uint16_t extend_idx;

			InToRead_t & in_to_read = intoread[id];

			read_idx = in_to_read.read_idx;
			seed_idx = in_to_read.seed_idx;
			extend_idx = in_to_read.extend_idx;


			
			Read_t & read_tmp = reads[read_idx];
			Seed_t & seed_tmp = read_tmp.seedInfo[seed_idx];

			// check if reverse complement
			if (seed_tmp.reverse){
				memcpy(in[i][j].read, read_tmp.pckSeq_rev, CEIL(read_tmp.readLen * BP_BIT, 8));
				memcpy(in[i][j].qs, read_tmp.pckQS_rev, CEIL(read_tmp.readLen * QS_BIT, 8));

				//cout<<"reversed: "<<std::bitset<8>(read_tmp.pckSeq_rev[0])<<" "<<std::bitset<8>(in[i][j].read[0])<<"\n";
				//cout<<"reversed QS: "<<std::bitset<8>(read_tmp.pckQS_rev[0])<<" "<<std::bitset<8>(in[i][j].qs[0])<<"\n";
			}else{
				memcpy(in[i][j].read, read_tmp.pckSeq, CEIL(read_tmp.readLen * BP_BIT, 8));
				memcpy(in[i][j].qs, read_tmp.pckQS, CEIL(read_tmp.readLen * QS_BIT, 8));

				//cout<<"normal: "<<std::bitset<8>(read_tmp.pckSeq[0])<<" "<<std::bitset<8>(in[i][j].read[0])<<"\n";
				//cout<<"normal QS: "<<std::bitset<8>(read_tmp.pckQS[0])<<" "<<std::bitset<8>(in[i][j].qs[0])<<"\n";
			}

			

			memcpy(in[i][j].ref, seed_tmp.pckRef[extend_idx], CEIL(REF_DEF_LEN * BP_BIT, 8));

			in[i][j].readLen = read_tmp.readLen;
			in[i][j].refLen = REF_DEF_LEN;

			in[i][j].isPad[0] = 0;


		}

		for (uint32_t j = 0; j < LATENCY*2; j++) {
			if (j < LATENCY*2 - 1){
				in[i][part_size[i]+j].id = rds_offset + part_size[i] + j;
				in[i][part_size[i]+j].isPad[0] = 1;

				in[i][part_size[i]+j].readLen = 101;

			}else{
				in[i][part_size[i]+j].id = rds_offset + part_size[i] + j;
				in[i][part_size[i]+j].isPad[0] = 2;

				in[i][part_size[i]+j].readLen = 101;
			}
		}

		rds_offset += part_size[i];
	}
	printf("OK");


		
	// align results (FPGA)
	//gettimeofday(&tv1, NULL);
	for (uint8_t i = 0; i < N_KRNL; i++){
		out[i] = new out_t [part_size[i]];
	    if (!out[i]) {
	    	 fprintf(stderr, "error: unable to allocate memory!\n");
	    	 exit(1);
	    }
	    memset(out[i], 0, part_size[i] * sizeof(out_t));
	}

	actionAlign(in,
				part_size,
				out,
				fpga_var.engine,
				bmk);



	printf("\t parsing results ... "); fflush(stdout);
	gettimeofday(&tv1, NULL);

	uint32_t total_aligned = 0;
	parse_result(reads, intoread, part_size, out, total_aligned);

	gettimeofday(&tv2, NULL);
	
	bmk->aligned_cnt = bmk->aligned_cnt + total_aligned;

	printf("[%.2f s] OK, Total Aligned: %d, Total UnAligned %d total %d \n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
		(double) (tv2.tv_sec - tv1.tv_sec),
		total_aligned,
		reads.size() - total_aligned,
		reads.size());

	// cleanup
	for (int i = 0; i < N_KRNL; i++) {
		delete[] in[i];
		delete[] out[i];
	}

}


void actionAlign(in_t **in,
					uint32_t *part_size,
					out_t **out,
					max_engine_t *engine, 
					bmk_t * bmk){

	struct timeval  tv1, tv2;
	double time_taken;


	uint64_t n_bytes_in[N_KRNL];
	uint64_t n_bytes_out[N_KRNL];

	for (uint8_t i = 0; i < N_KRNL; i++) {
		n_bytes_in[i] = (part_size[i] + LATENCY*2) * sizeof(in_t);
		n_bytes_out[i] = part_size[i] * sizeof(out_t);
	}

//cout<<"out size: "<<(int)(1 * sizeof(out_t))<<"\n";

	SmithAlpha_Align_actions_t * smith_align = new SmithAlpha_Align_actions_t;

	smith_align->param_interleave = LATENCY;
	smith_align->param_bramDepth = LATENCY * (REF_DEF_LEN + READ_DEF_LEN-1);

	smith_align->param_thresholdScoreIn = (int)(20 + 8.0 * log(HARD_CODE_READ_LEN));

	smith_align->param_MA = (uint8_t) 2;
	smith_align->param_NP = (uint8_t) 1;

	smith_align->param_MX = (uint8_t) 6;
	smith_align->param_MN = (uint8_t) 2;

	smith_align->param_RFG_OPEN = (uint8_t) 5;//should be 8
	smith_align->param_RFG_EXTEND = (uint8_t) 3;

	smith_align->param_RDG_OPEN = (uint8_t) 5;
	smith_align->param_RDG_EXTEND = (uint8_t) 3;

	smith_align->param_nBytesInput = n_bytes_in;
	smith_align->param_nBytesOutput = n_bytes_out;

#if N_KRNL > 0
	smith_align->instream_seqIn0 = (uint64_t*)in[0];
	smith_align->outstream_outRsltOut0 = (uint64_t*)out[0];
#endif

#if N_KRNL > 1
	smith_align->instream_seqIn1 = (uint64_t*)in[1];
	smith_align->outstream_outRsltOut1 = (uint64_t*)out[1];
#endif

#if N_KRNL > 2
	smith_align->instream_seqIn2 = (uint64_t*)in[2];
	smith_align->outstream_outRsltOut2 = (uint64_t*)out[2];
#endif

	gettimeofday(&tv1, NULL);
	SmithAlpha_Align_run(engine, smith_align);
    gettimeofday(&tv2, NULL);


    time_taken = ((double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
          (double) (tv2.tv_sec - tv1.tv_sec));

    printf("**FPGA TIME:** [%.2f s]\n", time_taken);
    fflush(stdout);

	bmk->process_time = bmk->process_time + time_taken;
	
	delete smith_align;
}


#endif /* ALIGN_H_ */