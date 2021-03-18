#include "format.hpp"
#include "align.hpp"

#include "fm-sa-def.h"


int getChrsIdx(uint32_t pos, std::vector<chr_t> &chrs){
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


void writePosDirectly(FILE *fp_di, 
                        std::vector<Read_t> &reads,
                        std::vector<InToRead_t> &intoread,
                        std::vector<chr_t> &chrs,
                        char * buffer){
                        
	uint64_t bytes = 0;


    uint32_t tmp_val = 0;
    int chr_idx;

    for (uint32_t i = 0; i < reads.size(); i++) {		

		if ((bytes + 1024) > BUFF_SIZE) {
            writeFile(fp_di, buffer, bytes);
            bytes = 0;
        }

        Read_t &read = reads[i];

        if (read.q.size() > 0){
            writeAlignedStr(bytes, buffer, (void *)(read.readId.c_str()), read.readId.size());
			buffer[bytes++] = '\n';

            uint32_t aligned_num = read.q.size();

            for (uint32_t h = 0; h < aligned_num; h++){
                AlignedInfo_t tmp = read.q.top();

                uint16_t seed_idx = tmp.seed_idx;
                uint16_t extend_idx = tmp.extend_idx;

                uint64_t refStart = read.seedInfo[seed_idx].refStart[extend_idx] -
                                        read.seedInfo[seed_idx].rdStart;

                chr_idx = getChrsIdx(refStart, chrs);
                writeAlignedStr(bytes, buffer, (void *)(chrs[chr_idx].name.c_str()+1), chrs[chr_idx].name.length()-2);
                buffer[bytes++] = ' ';

                tmp_val = refStart - chrs[chr_idx].begin;
                writeAlignedVal(bytes, buffer, &tmp_val, 10);
				buffer[bytes++] = ' ';

                read.q.pop();

            }
            buffer[bytes++] = '\n';
        }
    }

	if (bytes > 0) {
        writeFile(fp_di, buffer, bytes);
    }	
}