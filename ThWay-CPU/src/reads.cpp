/*
 * reads.cpp
 *
 *  Created on: Feb 3, 2019
 *      Author: hn915
 */

#include "reads.h"


void packSymbols(char *sym, uint8_t *pck, uint8_t len);
inline void setVal(uint8_t *pck, uint32_t idx, uint8_t val);

// get remaining entry symbols if '@' is a quality score
uint64_t getOneMoreEntry(FILE *fp, char *buffer, uint64_t len)
{
    long int f_pos = ftell(fp);
    int cnt = 0;
    int c;

    // test if entry
    while ((c = fgetc(fp)) != EOF) {
        cnt += 1;
        if (c == '\n') {
            c = fgetc(fp);

            // It can be some random stuff?
            if (c != '@' && c != EOF)
                cnt = 0;
            break;
        }
    }

    // reset file pointer
    fseek(fp, f_pos, SEEK_SET);

    // read remaining positions
    for (int i = 0; i < cnt; i++)
        buffer[len++] = fgetc(fp);

    return len;
}




void loadReads(FILE *fp, std::vector<read2Bit_t> &reads2Bit, 
		char *buffer, uint64_t size_r, uint64_t *bytes, bool r_ctrl){

    struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);

	reads2Bit.clear();

	if (r_ctrl == true && size_r > 0) {

		// read a bunch of things from file first
		if (fread(buffer, size_r, 1, fp) != 1) {
			fprintf(stderr, "error: unable to read file!\n");
			exit(1);
		}

		// Then read until next record starts or eof
		int c;
		while ((c = fgetc(fp)) != EOF) {
			if (c == '@') {
				ungetc(c, fp);
				size_r = getOneMoreEntry(fp, buffer, size_r);
				break;
			}
			buffer[size_r++] = (char)c;
		}

		// update bytes read
		*bytes += size_r;


		// parse the buffer
		read2Bit_t tmp2Bit;

		//char tmp_seq[MAX_READ_LEN+1];
		std::string tmp_at_line;
		std::string tmp_q_score;
		uint32_t tmp_seq_len;

		tmp2Bit.is_align = false;
		tmp2Bit.n_seedhits = 0;


		uint64_t i = 0;
		uint64_t s_position;
		uint64_t s_position_seq;
		int s_len;

		while (i < size_r) {
			bool has_N = false;
			//tmp2Bit.has_N = false;
			//tmp3Bit.has_N = true;

			//// Step 1: get meta data line, which is <@r337>
			s_position = i;
			s_len = 0;

			while (buffer[i++] != '\n') {
				s_len += 1;
			}

			tmp_at_line.assign(&buffer[s_position], s_len);


			////Step 2: get the sequence
			//bool unknown = false;
			s_position = i;
			s_position_seq = i;
			s_len = 0;

			while (buffer[i] != '\n') {
				i++;
				//if (buffer[i] == 'N')
					//has_N = true;
				s_len += 1;
			}

			i += 1;

			//check if the length of the sequence is longer than the allocated size (not likely to happen)
			s_len = s_len > MAX_READ_LEN ? MAX_READ_LEN : s_len;

			//do the str copying much later
			//		memcpy(tmp_seq, &buffer[s_position], s_len * sizeof(char));
			//		tmp_seq[s_len] = '\0';
			tmp_seq_len = s_len;


			////Step 3: get the line with '+' marker
			//tmp.plus_data = buffer[i];
			while (buffer[i++] != '\n');

			////Step 4: get the quality score
			s_position = i;
			s_len = 0;

			while (buffer[i++] != '\n') {
				s_len += 1;
			}

			tmp_q_score.assign(&buffer[s_position], s_len);

			if (has_N == false){

				tmp2Bit.at_line = tmp_at_line;
				//tmp2Bit.seq = tmp_seq;
				tmp2Bit.seq_len = tmp_seq_len;
				tmp2Bit.q_score = tmp_q_score;

				memcpy(tmp2Bit.seq, &buffer[s_position_seq], tmp_seq_len * sizeof(char));
				tmp2Bit.seq[tmp_seq_len] = '\0';

				reads2Bit.push_back(tmp2Bit);

			}
		}

		int n_threads = omp_get_max_threads();
		//int n_threads = 7;

		// change the symbol to binary
#pragma omp parallel for num_threads(n_threads)
		for (uint32_t i = 0; i < reads2Bit.size(); i++) {
			memset(reads2Bit[i].pck_sym, 0, CEIL(reads2Bit[i].seq_len, FM_BP_RANGE) * sizeof(uint8_t));
			packSymbols(reads2Bit[i].seq, reads2Bit[i].pck_sym, reads2Bit[i].seq_len);

			//pack the sub symbol
			uint8_t divlen =  (reads2Bit[i].seq_len) / MISMATCH_DIV;
			//cout<<"\n\ncur_seq_name: "<<(int )(reads2Bit[i].seq_len) <<"\n";

			for (uint8_t n = 0; n < MISMATCH_DIV;  n++){
				uint8_t offset = divlen * n;
				uint8_t curlen = (n == MISMATCH_DIV - 1 )? 
									reads2Bit[i].seq_len - divlen * (MISMATCH_DIV-1):
									divlen;
				reads2Bit[i].seq_sub_len[n] = curlen;
				reads2Bit[i].seq_sub_pos[n] = offset;
				
				memset(reads2Bit[i].pck_sub_sym[n], 0, (CEIL(MAX_READ_LEN, FM_BP_RANGE)/4 + 1) * sizeof(uint8_t));
				//cout<<"LEN: "<<(int)curlen<<"\n";
				packSymbols(reads2Bit[i].seq + offset, reads2Bit[i].pck_sub_sym[n], curlen);
				//cout<<"\n";
			}
		}
	}

	gettimeofday(&tv2, NULL);
    printf("\nFINISH ---> Loading reads [%.2f s]\n\n",(double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
			(double) (tv2.tv_sec - tv1.tv_sec));

	return;
}

// pack symbols
void packSymbols(char *sym, uint8_t *pck, uint8_t len){
  for (uint8_t i = 0; i < len; i++) {
	// backward search, so we start from the end
    switch(sym[len-1-i]) {
    case 'A': setVal(pck, i, 0); break;
    case 'C': setVal(pck, i, 1); break;
    case 'G': setVal(pck, i, 2); break;
    case 'T': setVal(pck, i, 3); break;
    case 'a': setVal(pck, i, 0); break;
    case 'c': setVal(pck, i, 1); break;
    case 'g': setVal(pck, i, 2); break;
    case 't': setVal(pck, i, 3); break;
    default : setVal(pck, i, 0);
    }
  }
}

inline void setVal(uint8_t *pck, uint32_t idx, uint8_t val){
  uint8_t tmp = val << ((idx * FM_BP_BIT) % (sizeof(uint8_t)*8));
  pck[idx/FM_BP_RANGE] |= tmp;
}

uint64_t writeReadsN(FILE *fp, std::vector<read3Bit_t> &reads, char *buffer){
    uint64_t bytes = 0;

	for (uint32_t i = 0; i < reads.size(); i++) {
		// write buffer before overflow
		if ((bytes + 512) > BUFF_SIZE) {
			writeFile(fp, buffer, bytes);
			bytes = 0;
		}

		// write meta data
		memcpy(&buffer[bytes], reads[i].at_line.c_str(), reads[i].at_line.size());
		bytes += reads[i].at_line.size();
		buffer[bytes++] = '\n';

		// write seqeunce data
		memcpy(&buffer[bytes], reads[i].seq, reads[i].seq_len);
		bytes += reads[i].seq_len;
		buffer[bytes++] = '\n';

		// write strand
		buffer[bytes++] = '+';
		buffer[bytes++] = '\n';

		// write quality scores
		memcpy(&buffer[bytes], reads[i].q_score.c_str(), reads[i].q_score.size());
		bytes += reads[i].q_score.size();
		buffer[bytes++] = '\n';
	}
	// write the rest if there are any
	if (bytes > 0) {
		writeFile(fp, buffer, bytes);
		bytes = 0;
	}

	return reads.size();
}


uint64_t writeReads(FILE *fp, std::vector<read2Bit_t> &reads, char *buffer){
    uint64_t bytes = 0;
    uint64_t aligned_cnt = 0;
    uint64_t unaligned_cnt = 0;

    // write hits to buffer
    for (uint32_t i = 0; i < reads.size(); i++) {
    	if (reads[i].n_seedhits == 0) {
			// write buffer before overflow
			if ((bytes + 512) > BUFF_SIZE) {
				writeFile(fp, buffer, bytes);
				bytes = 0;
			}

			// write meta data
			memcpy(&buffer[bytes], reads[i].at_line.c_str(), reads[i].at_line.size());
			bytes += reads[i].at_line.size();
			buffer[bytes++] = '\n';

			// write seqeunce data
			memcpy(&buffer[bytes], reads[i].seq, reads[i].seq_len);
			bytes += reads[i].seq_len;
			buffer[bytes++] = '\n';

			// write strand
			buffer[bytes++] = '+';
			buffer[bytes++] = '\n';

			// write quality scores
			memcpy(&buffer[bytes], reads[i].q_score.c_str(), reads[i].q_score.size());
			bytes += reads[i].q_score.size();
			buffer[bytes++] = '\n';

			unaligned_cnt = unaligned_cnt + 1;

    	}
    	else{
            aligned_cnt = aligned_cnt + 1;

			// write as SAM file
			
    	}
    }


    if (bytes > 0) {
        writeFile(fp, buffer, bytes);
    }

    return aligned_cnt;
}
