/*
 * reads.cpp
 *
 *  Created on: Feb 3, 2019
 *      Author: hn915
 */

#include "reads.h"
#include <bitset>


void packSymbols(char *sym, uint8_t *pck, uint8_t len);
inline void setVal(uint8_t *pck, uint32_t idx, uint8_t val);

void packSymQS(char *sym, const char *qs, uint8_t * pck, uint8_t len);
inline void setSymQSVal(uint8_t *pck, uint32_t idx, uint8_t symVal, char qsVal);

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




void loadReads(FILE *fp, std::vector<read2Bit_t> &reads2Bit, //std::vector<read3Bit_t> &reads3Bit,
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
		//read3Bit_t tmp3Bit;

		//char tmp_seq[MAX_READ_LEN+1];
		std::string tmp_at_line;
		std::string tmp_q_score;
		uint32_t tmp_seq_len;

		tmp2Bit.is_align = false;
		tmp2Bit.n_hits = 0;
		//tmp2Bit.low = 0;
		//tmp2Bit.high = 0;

		//tmp3Bit.is_align = false;
		//tmp3Bit.low = 0;
		//tmp3Bit.high = 0;

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
				if (buffer[i++] == 'N')
					has_N = true;
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

			}/*else{

				tmp3Bit.at_line = tmp_at_line;
				//tmp3Bit.seq = tmp_seq;
				tmp3Bit.seq_len = tmp_seq_len;
				tmp3Bit.q_score = tmp_q_score;

				memcpy(tmp3Bit.seq, &buffer[s_position_seq], tmp_seq_len * sizeof(char));
				tmp3Bit.seq[tmp_seq_len] = '\0';

				reads3Bit.push_back(tmp3Bit);

			}*/

		}

		int n_threads = omp_get_max_threads();
		//int n_threads = 7;

		// change the symbol to binary
		/*#pragma omp parallel for num_threads(n_threads)
		for (uint32_t i = 0; i < reads2Bit.size(); i++) {
			memset(reads2Bit[i].pck_sym, 0, CEIL(reads2Bit[i].seq_len, FM_BP_RANGE) * sizeof(uint8_t));
			packSymbols(reads2Bit[i].seq, reads2Bit[i].pck_sym, reads2Bit[i].seq_len);
		}*/

		// change the symbol and quality score to binary
#pragma omp parallel for num_threads(n_threads)
		for (uint32_t i = 0; i < reads2Bit.size(); i++) {
			memset(reads2Bit[i].pck_sym_qs, 0, reads2Bit[i].seq_len * sizeof(uint8_t));
			packSymQS(reads2Bit[i].seq, reads2Bit[i].q_score.c_str(), reads2Bit[i].pck_sym_qs, reads2Bit[i].seq_len);
		}

		//change the quality score to binary
	}
	gettimeofday(&tv2, NULL);
	printf("OK Load Reads Time [%.2f s]\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
				(double) (tv2.tv_sec - tv1.tv_sec));
	return;
}

void packSymQS(char *sym, const char *qs, uint8_t * pck, uint8_t len){
	for (uint8_t i = 0; i < len; i++) {
	// backward search, so we start from the end
		switch(sym[len-1-i]) {
		case 'A': setSymQSVal(pck, i, 0, qs[len-1-i]); break;
	    case 'C': setSymQSVal(pck, i, 1, qs[len-1-i]); break;
	    case 'G': setSymQSVal(pck, i, 2, qs[len-1-i]); break;
	    case 'T': setSymQSVal(pck, i, 3, qs[len-1-i]); break;
		case 'a': setSymQSVal(pck, i, 0, qs[len-1-i]); break;
	    case 'c': setSymQSVal(pck, i, 1, qs[len-1-i]); break;
	    case 'g': setSymQSVal(pck, i, 2, qs[len-1-i]); break;
	    case 't': setSymQSVal(pck, i, 3, qs[len-1-i]); break;
	    default : setVal(pck, i, 0);
		}
	}
}

inline void setSymQSVal(uint8_t *pck, uint32_t idx, uint8_t symVal, char qsVal){
	uint8_t tmp = ((uint8_t) qsVal - 33) << 2;
	pck[idx] = symVal | tmp;
	//std::cout << "======= " << std::bitset<8>((char)pck[idx]) << qsVal<<"\n";
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

/*void writeReadsN(FILE *fp, std::vector<read3Bit_t> &reads, char *buffer){
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
		buffer[bytes++] = '\n'inline void writeAlignedStr(uint64_t &bytes, char *buffer, void *a, uint64_t a_bytes){ 
    
    char * a_char = (char *)a;
    uint16_t end_idx = a_bytes;
    
    // custom loop to get the first space
    for (uint64_t i = 0; i < a_bytes; i++){
        if (a_char[i] == ' '){
            end_idx = i;
            break;
        }
    }
    memcpy(&buffer[bytes], a, sizeof(char)*end_idx);
    bytes += end_idx;
}

}*/


uint64_t writeReads(FILE *fp, std::vector<read2Bit_t> &reads, char *buffer){
    uint64_t bytes = 0;
    uint64_t aligned_cnt = 0;

    // write hits to buffer
    for (uint32_t i = 0; i < reads.size(); i++) {
    	if (reads[i].hits.size() == 0 ) {
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

    	}else{
            aligned_cnt = aligned_cnt + 1;
    	}
    }


    if (bytes > 0) {
        writeFile(fp, buffer, bytes);
    }

    return aligned_cnt;
}