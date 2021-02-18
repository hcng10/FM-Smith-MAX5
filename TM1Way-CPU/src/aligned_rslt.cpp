#include "aligned_rslt.h"

inline void writeAligned(uint64_t &bytes, char *buffer, void *a, uint64_t a_bytes){   

    memcpy(&buffer[bytes], a, a_bytes);
	bytes += a_bytes;
}

inline void writeAlignedVal(uint64_t &bytes, char *buffer, void *a, uint64_t a_bytes){

    char tmp_buffer [a_bytes];
    memset(tmp_buffer, 0, sizeof(tmp_buffer));

    // convert to string 
    sprintf(tmp_buffer, "%u", *(uint32_t*)a);
    
    // write each char string element
    writeAligned(bytes, buffer, (void *)tmp_buffer, sizeof(char));
    for (uint8_t i = 1; i<sizeof(tmp_buffer); i++){
        if (tmp_buffer[i] != 0){
            writeAligned(bytes, buffer,  (void *) (tmp_buffer + i), sizeof(char));
        }
    }
}

inline void writeAlignedStr(uint64_t &bytes, char *buffer, void *a, uint64_t a_bytes){ 
    
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

void writeAligned(FILE *fp, std::vector<read2Bit_t> &reads, char *buffer){

	uint64_t bytes = 0;

	uint32_t tmp_val = 0;

	// write hits to buffer
    for (uint32_t i = 0; i < reads.size(); i++) {
		if (reads[i].n_hits > 0 ) {

			read2Bit_t read = reads[i];
			uint8_t n_hits = reads[i].n_hits;

			// write buffer before overflow
			if ((bytes + 512) > BUFF_SIZE) {
				writeFile(fp, buffer, bytes);
				bytes = 0;
			}
			// FORMAT: ReadID NumberofHit rev_cmpt? low0 high0 mispos(I) missym(I) qs(I) mispos(II) missym(II) qs(II) rev_cmpt? ...
			// write read name
			writeAlignedStr(bytes, buffer, (void *)(read.at_line.c_str()+1), read.at_line.size()-1);
			buffer[bytes++] = ' ';

			writeAlignedStr(bytes, buffer, (void *)(read.seq), read.seq_len);
			buffer[bytes++] = ' ';

			writeAlignedStr(bytes, buffer, (void *)(read.q_score.c_str()), read.seq_len);
			buffer[bytes++] = ' ';

			tmp_val = (uint32_t) n_hits;
			writeAlignedVal(bytes, buffer, &tmp_val, MATCH_BIT);
			buffer[bytes++] = ' ';

			// get all the hit details <32, one loop is OK already
			for (uint8_t h = 0; h < n_hits; h++){
				tmp_val = (read.hits[0].is_rev_cmpt >> (N_HITS-1-h)) & 1;
				writeAlignedVal(bytes, buffer, &tmp_val, 1);
				buffer[bytes++] = ' ';

				tmp_val = read.hits[0].low_sorted[N_HITS-1-h];
				writeAlignedVal(bytes, buffer, &tmp_val, POS_BIT);
				buffer[bytes++] = ' ';

				tmp_val = read.hits[0].high_sorted[N_HITS-1-h];
				writeAlignedVal(bytes, buffer, &tmp_val, POS_BIT);
				buffer[bytes++] = ' ';

				tmp_val = (uint32_t) read.hits[0].mis_pos1[N_HITS-1-h];
				writeAlignedVal(bytes, buffer, &tmp_val, POS_BIT);
				buffer[bytes++] = ' ';

				buffer[bytes++] = read.hits[0].mis_sym1[N_HITS-1-h];
				buffer[bytes++] = ' ';

				buffer[bytes++] = read.hits[0].qs1[N_HITS-1-h];
				buffer[bytes++] = ' ';	

				tmp_val = (uint32_t) read.hits[0].mis_pos2[N_HITS-1-h];
				writeAlignedVal(bytes, buffer, &tmp_val, POS_BIT);
				buffer[bytes++] = ' ';

				buffer[bytes++] = read.hits[0].mis_sym2[N_HITS-1-h];
				buffer[bytes++] = ' ';

				buffer[bytes++] = read.hits[0].qs2[N_HITS-1-h];
				buffer[bytes++] = ' ';		
				
			}
			buffer[bytes++] = '\n';
		}
	}

	if (bytes > 0) {
        writeFile(fp, buffer, bytes);
    }

	return;
}

inline int getChrsIdx(uint32_t pos, std::vector<chr_t> &chrs){
	if (chrs.size() == 0) return 0;

	uint32_t last = chrs.size() - 1;
    // Target index bigger or smaller than the indices 
    if (pos < chrs[0].begin) return 0;
	if (pos > chrs[last].begin) return last;

	uint16_t start = 0;
    uint16_t end = chrs.size();
	uint16_t mid = 0;


    while (start < end){
        mid = (start + end) / 2;

		if (pos == chrs[mid].begin){
			return mid;
		}

		if (pos < chrs[mid].begin){
			if (mid > 0 && pos > chrs[mid-1].begin){
				return mid-1;
			}
			end = mid;
		}
		else{
			if (mid < last && pos <  chrs[mid+1].begin){
				return mid;
			}
			start = mid + 1;
		}
	}
	return mid;
}


void writePosDirectly(FILE *fp_di, 
                        std::vector<read2Bit_t> &reads,
                        uint32_t *sai,
                        uint32_t *sai_rev,
                        std::vector<chr_t> &chrs,
                        char * buffer){
						
	uint64_t bytes = 0;

	uint32_t low_tmp = 0;
	uint32_t high_tmp = 0;
	uint8_t aligned_bck = 0;


    uint32_t sa_tmp = 0;
    uint32_t tmp_val = 0;
	
    //char tmp_char = '0';
    int chr_idx;
	int james_cnt = 0;

	for (uint32_t i = 0; i < reads.size(); i++) {		

		if ((bytes + BUFF_SIZE/3) > BUFF_SIZE) {
            writeFile(fp_di, buffer, bytes);
            bytes = 0;
        }
	

		read2Bit_t read = reads[i];

		if (read.n_hits > 0 ) {

			writeAlignedStr(bytes, buffer, (void *)(read.at_line.c_str()), read.at_line.size());
			buffer[bytes++] = '\n';
#if IS_JAMES == 1
			james_cnt = 0;
#endif
			

			for (uint8_t h = 0; h < read.hits.size(); h++){
				for (uint8_t hi = 0; hi < read.hits[h].n_hits; hi++){

					low_tmp = read.hits[h].low_sorted[N_HITS-1-hi];
					high_tmp = read.hits[h].high_sorted[N_HITS-1-hi];

					aligned_bck = (read.hits[h].is_aligned_bck >> (N_HITS-1-hi)) & 1;

					/*if (i >= 15659745){
						cerr<<"number of hit: "<<(int)read.hits[h].n_hits<<"\n";
						cerr<<high_tmp<<" "<<low_tmp<<" "<<(high_tmp - low_tmp)<<" "<<bytes<<" Buff size inside: "<<BUFF_SIZE<<"\n";
					}*/
					

					for (uint32_t j = 0; j < high_tmp - low_tmp && j < 144; j++){
						if (aligned_bck == 1){
							sa_tmp = sai[low_tmp + j] + 1;
						}else{
							sa_tmp = sai_rev[low_tmp + j] + 1 - (read.seq_len -1);
						}
						
						chr_idx = getChrsIdx(sa_tmp, chrs);
                    	writeAlignedStr(bytes, buffer, (void *)(chrs[chr_idx].name.c_str()+1), chrs[chr_idx].name.length()-2);
                   		buffer[bytes++] = ' ';

						tmp_val = sa_tmp - chrs[chr_idx].begin;
                    	writeAlignedVal(bytes, buffer, &tmp_val, POS_BIT);
						buffer[bytes++] = ' ';

					}
#if IS_JAMES == 1					
					james_cnt = james_cnt + 1;
					if (james_cnt > 2){
						james_cnt = 0;
						break;
					}
#endif
				}
			}
			buffer[bytes++] = '\n';

		}
	}
	
	if (bytes > 0) {
        writeFile(fp_di, buffer, bytes);
    }	
}