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
			// FORMAT: ReadID NumberofHit low0 high0 mispos0 missym0 qs0 low1 ...
			// write read name
			writeAlignedStr(bytes, buffer, (void *)(read.at_line.c_str()+1), read.at_line.size()-1);
			buffer[bytes++] = ' ';

			tmp_val = (uint32_t) n_hits;
			writeAlignedVal(bytes, buffer, &tmp_val, MATCH_BIT);
			buffer[bytes++] = ' ';

			// get all the hit details <32
			for (uint8_t h = 0; h < read.hits.size(); h++){
				for (uint8_t hi = 0; hi < read.hits[h].n_hits; hi++){

					tmp_val = read.hits[h].low_sorted[N_HITS-1-hi];
					writeAlignedVal(bytes, buffer, &tmp_val, POS_BIT);
					buffer[bytes++] = ' ';

					tmp_val = read.hits[h].high_sorted[N_HITS-1-hi];
					writeAlignedVal(bytes, buffer, &tmp_val, POS_BIT);
					buffer[bytes++] = ' ';

					tmp_val = (uint32_t) read.hits[h].mis_pos_sorted[N_HITS-1-hi];
					writeAlignedVal(bytes, buffer, &tmp_val, POS_BIT);
					buffer[bytes++] = ' ';

					buffer[bytes++] = read.hits[h].mis_sym_sorted[N_HITS-1-hi];
					buffer[bytes++] = ' ';

					buffer[bytes++] = read.hits[h].qs_sorted[N_HITS-1-hi];
					buffer[bytes++] = ' ';				
				}
			}
			buffer[bytes++] = '\n';
		}
	}

	if (bytes > 0) {
        writeFile(fp, buffer, bytes);
    }

	return;
}