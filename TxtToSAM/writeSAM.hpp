#ifndef WRITESAM_H
#define WRITESAM_H

#include <stdint.h>
#include <omp.h>
#include <type_traits>


#include "format.h"
#include "fm-sa-def.h"
#include "file_op.h"

#define POS_BIT 10
#define QS_BIT 3
#define MATCH_BIT 5
#define ED_DIST_BIT 2

#define ALIGN_REPORT_NUM 16


inline void writeSAM(uint64_t &bytes, char *buffer, void *a, uint64_t a_bytes){   

    memcpy(&buffer[bytes], a, a_bytes);
	bytes += a_bytes;
}

inline void writeSAMCStr(uint64_t &bytes, char *buffer, void *a, uint64_t a_bytes){ 
    
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


inline void writeSAMTab(uint64_t &bytes, char *buffer){
    
    char t = '\t';

    memcpy(&buffer[bytes], &t, 1);
	bytes += 1;
}


inline void writeSAMVal(uint64_t &bytes, char *buffer, void *a, uint64_t a_bytes){

    char tmp_buffer [a_bytes];
    memset(tmp_buffer, 0, sizeof(tmp_buffer));

    // convert to string 
    sprintf(tmp_buffer, "%u", *(uint32_t*)a);
    
    // write each char string element
    writeSAM(bytes, buffer, (void *)tmp_buffer, sizeof(char));
    for (uint8_t i = 1; i<sizeof(tmp_buffer); i++){
        if (tmp_buffer[i] != 0){
            writeSAM(bytes, buffer,  (void *) (tmp_buffer + i), sizeof(char));
        }
    }
}

template <typename hit_t>
inline void setSAMOptional(uint64_t &bytes, 
                            char *buffer, 
                            read2Bit_t & read, 
                            hit_t * hit,
                            uint32_t & suppress_aligned,
                            bool is_om){

    // MD:Z:<S> For aligned reads, 
    //          <S> is a string representation of the 
    //          mismatched reference bases in the alignment
    // XA:i:<N> Aligned read belongs to stratum <N>.
    // XM:i:<N> the read’s alignments were suppressed
    // NM:i:<V> edit distance
    
    char tmp_char = '0';
    char tmp_cstr[] = "MD:Z:";
    uint32_t tmp_val = read.seq_len;

    writeSAM(bytes, buffer, &tmp_cstr, sizeof(tmp_cstr)-1);
    if (is_om == 1){
        hit_om_t * hit_om = (hit_om_t*) hit;
        tmp_val = hit_om->mis_pos;
        writeSAMVal(bytes, buffer, &tmp_val, POS_BIT);

        tmp_char = hit_om->mis_sym;
        writeSAM(bytes, buffer, &tmp_char, 1);

        if (read.seq_len > tmp_val){
            tmp_val = read.seq_len - tmp_val -1;
            writeSAMVal(bytes, buffer, &tmp_val, MATCH_BIT);
        }
    }else{
        hit_tm_t * hit_tm = (hit_tm_t*) hit;
        //check which mispos is smaller
        uint8_t small_pos = 0;
        uint8_t big_pos = 0;

        char small_sym;
        char big_sym;

        if (hit_tm->mis_pos0 < hit_tm->mis_pos1){
            small_pos = hit_tm->mis_pos0;
            big_pos = hit_tm->mis_pos1;

            small_sym = hit_tm->mis_sym0;
            big_sym = hit_tm->mis_sym1;
        }else{
            small_pos = hit_tm->mis_pos1;
            big_pos = hit_tm->mis_pos0;

            small_sym = hit_tm->mis_sym1;
            big_sym = hit_tm->mis_sym0;
        }


        tmp_val = small_pos;
        writeSAMVal(bytes, buffer, &tmp_val, POS_BIT);

        tmp_char = small_sym;
        writeSAM(bytes, buffer, &tmp_char, 1);

        if (small_pos + 1 != big_pos){
            tmp_val = big_pos - small_pos - 1;
            writeSAMVal(bytes, buffer, &tmp_val, POS_BIT);
        }

        tmp_char = big_sym;
        writeSAM(bytes, buffer, &tmp_char, 1);

        if (read.seq_len > tmp_val){
            tmp_val = read.seq_len - big_pos -1;
            writeSAMVal(bytes, buffer, &tmp_val, MATCH_BIT);
        }

    }
    writeSAMTab(bytes, buffer);

    tmp_char = '0';

    strcpy(tmp_cstr, "XA:i:");
    writeSAM(bytes, buffer, &tmp_cstr, sizeof(tmp_cstr)-1);
    writeSAM(bytes, buffer, &tmp_char, 1);
    writeSAMTab(bytes, buffer);


    if (suppress_aligned > 0){    
        strcpy(tmp_cstr, "XM:i:");
        writeSAMVal(bytes, buffer, &suppress_aligned, POS_BIT);
        writeSAMTab(bytes, buffer);
    }

    
    strcpy(tmp_cstr, "NM:i:");
    writeSAM(bytes, buffer, &tmp_cstr, sizeof(tmp_cstr)-1);
    if (is_om == 1){
        tmp_char = '1';
    }
    writeSAM(bytes, buffer, &tmp_char, 1);
}


template <typename hit_t>
inline void setSAMFlag(uint64_t &bytes, char *buffer, hit_t & hit){

    int flag = 0;

    // 0x10
    if (hit.is_rev_cmpt == 1){
        flag += 16;
    }

    writeSAMVal(bytes, buffer, &flag, 4);
}

inline char getReverseComplement(char b){

    char rtn_char;
    switch(b){
        case 'A':
            rtn_char = 'T';
            break;
        case 'C':
            rtn_char = 'G';
            break;
        case 'G':
            rtn_char = 'C';
            break;
        case 'T':
            rtn_char = 'A';
            break;
        
        case 'a':
            rtn_char = 't';
            break;
        case 'c':
            rtn_char = 'g';
            break;
        case 'g':
            rtn_char = 'c';
            break;
        case 't':
            rtn_char = 'a';
            break;
        default:
            rtn_char = 'N';
    }
    return rtn_char;
}

void getReverseStr(string& str) 
{ 
    int n = str.length(); 
  
    // Swap character starting from two 
    // corners 
    for (int i = 0; i < n / 2; i++) 
        swap(str[i], str[n - i - 1]); 
} 

int getChrsIdx(uint32_t pos, std::vector<chr_t> &chrs){

    if (chrs.size() == 0) return 0;
    // binary search
    uint16_t chrs_num = chrs.size();
    uint16_t mid = (chrs_num  - 1)/ 2;
    
    uint16_t start = 0;
    uint16_t end = chrs_num - 1;

    uint16_t sel = 0;

    while(start <= end){

        if (pos > chrs[mid].begin){
            sel = mid;

            start = mid + 1;
            end = end;
            mid = (end - start) / 2 + start; 
        }
        else if (pos < chrs[mid].begin){
            start = start;
            end = (mid == 0)? mid: mid - 1;
            mid = (end - start) / 2 + start; 

            sel = end;
        }
        else if (pos == chrs[mid].begin){
            return mid;
        }

        if (start >= end){
            if (pos >= chrs[end].begin){
                return end;
            }
            else{
                return sel;
            }
        }
    }
    return -1;
}

template <typename hit_t>
void writeSAMFields(uint64_t &bytes, 
                    char *   buffer,
                    std::vector<chr_t> &chrs,
                    read2Bit_t read,
                    hit_t * hit,
                    uint32_t & suppress_aligned,
                    uint32_t * aligned_sai_b,
                    uint32_t * aligned_sai_f,
                    bool is_b,
                    uint32_t sai_i,
                    bool is_om){
                    
    uint32_t tmp_val = 0;
    char tmp_char = '0';
    char seq_revc[MAX_READ_LEN+1];
    int chr_idx;
    //cerr<<"chk_pt1\n";

    // 1. QNAME
    writeSAMCStr(bytes, buffer, (void *)(read.at_line.c_str()), read.at_line.size());
    writeSAMTab(bytes, buffer);
    //cerr<<"chk_pt2\n";

    // 2. Write Flag
    if (is_om){
        hit_om_t * hit_om = (hit_om_t*) hit;
        setSAMFlag<hit_om_t>(bytes, buffer, * hit_om);
    }else{
        hit_tm_t * hit_tm = (hit_tm_t*) hit;
        setSAMFlag<hit_tm_t>(bytes, buffer, * hit_tm);
    }
    writeSAMTab(bytes, buffer);
    //cerr<<"chk_pt3\n";

    // 3. Search and Write RNAME
    if (chrs.size() == 0){
        tmp_char = '*';
        writeSAM(bytes, buffer, &tmp_char, 1);
    } else{
       if (is_b){
            chr_idx = getChrsIdx(aligned_sai_b[sai_i], chrs);
        }else{
            chr_idx = getChrsIdx(aligned_sai_f[sai_i]  - (read.seq_len - 1), chrs);
        }

        writeSAMCStr(bytes, buffer, (void *)(chrs[chr_idx].name.c_str()+1), chrs[chr_idx].name.length()-2);
    }
    writeSAMTab(bytes, buffer);   
    //cerr<<"chk_pt4\n";

    // 4. Write POS
    if (is_b){
        tmp_val = aligned_sai_b[sai_i] - chrs[chr_idx].begin;
        writeSAMVal(bytes, buffer, &tmp_val, POS_BIT);
    }
    else{
        tmp_val = aligned_sai_f[sai_i] - (read.seq_len - 1) + chrs[chr_idx].begin;
        writeSAMVal(bytes, buffer, &tmp_val, POS_BIT);
    }
    writeSAMTab(bytes, buffer);
    //cerr<<"chk_pt5\n";

    // 5. write QS
    if (is_om){
        hit_om_t * hit_om = (hit_om_t*) hit;
        tmp_val = (uint32_t) hit_om->qs - 33;
    }else{
        hit_tm_t * hit_tm = (hit_tm_t*) hit;
        tmp_val = (uint32_t) (hit_tm->qs0 - 33) + (hit_tm->qs1 - 33);
    }
    //tmp_val = (uint32_t) hit.qs - 33;
    writeSAMVal(bytes, buffer, &tmp_val, QS_BIT);
    writeSAMTab(bytes, buffer);
    //cerr<<"chk_pt6\n";
    
    // 6. write CIGAR
    tmp_val = read.seq_len;
    tmp_char = 'M';
    writeSAMVal(bytes, buffer, &tmp_val, MATCH_BIT);
    writeSAM(bytes, buffer, &tmp_char, 1);
    writeSAMTab(bytes, buffer);
    //cerr<<"chk_pt7\n";

    // 7. RNEXT (pair-end)
    tmp_char = '*';
    writeSAM(bytes, buffer, &tmp_char, 1);
    writeSAMTab(bytes, buffer);
    //cerr<<"chk_pt8\n";

    // 8. PNEXT (pair-end)
    tmp_char = '0';
    writeSAM(bytes, buffer, &tmp_char, 1);
    writeSAMTab(bytes, buffer);
    //cerr<<"chk_pt9\n";

    // 9. TLEN (pair-end)
    writeSAM(bytes, buffer, &tmp_char, 1);
    writeSAMTab(bytes, buffer);
    //cerr<<"chk_pt10\n";

    // 10. SEQ
    bool is_rev_cmpt;
    if (is_om){
        hit_om_t * hit_om = (hit_om_t*) hit;
        is_rev_cmpt = hit_om->is_rev_cmpt;
    }else{
        hit_tm_t * hit_tm = (hit_tm_t*) hit;
        is_rev_cmpt = hit_tm->is_rev_cmpt;
    }

    if (!is_rev_cmpt){
        writeSAMCStr(bytes, buffer, (void *)(read.seq.c_str()), read.seq_len);
    }else{
        for (uint32_t i = 0; i < read.seq_len; i++){
            tmp_char = read.seq[read.seq_len-1-i];
            seq_revc[i] = getReverseComplement(tmp_char);          
        }
        writeSAM(bytes, buffer, &seq_revc, read.seq_len);
    }
    writeSAMTab(bytes, buffer);            
    
    // 11. QUAL
    if (is_rev_cmpt){
        getReverseStr(read.q_score);
    }
    writeSAM(bytes, buffer, (void *)(read.q_score.c_str()), read.seq_len);
    writeSAMTab(bytes, buffer);


    // 12. Optional 
    if (is_om){
        hit_om_t * hit_om = (hit_om_t*) hit;
        setSAMOptional<hit_om_t>(bytes, buffer, read, hit_om, suppress_aligned, is_om);
    }else{
        hit_tm_t * hit_tm = (hit_tm_t*) hit;
        setSAMOptional<hit_tm_t>(bytes, buffer, read, hit_tm, suppress_aligned, is_om);
    }
    tmp_char = '\n';
    writeSAM(bytes, buffer, &tmp_char, 1);
}

void convertSAM(FILE *fp_sam, 
                char *buffer, 
                std::vector<read2Bit_t> &reads, 
                uint32_t *sai, 
                uint32_t * sai_rev, 
                std::vector<chr_t> &chrs, 
                bool is_om){

    uint64_t bytes = 0;
    uint32_t aligned_sai_b[ALIGN_REPORT_NUM];
    uint32_t aligned_sai_f[ALIGN_REPORT_NUM];

    uint32_t aligned_num_b = 0;
    uint32_t aligned_num_f = 0;
    uint32_t suppress_aligned = 0;

    for (uint32_t i = 0; i < reads.size(); i++) {
        read2Bit_t read = reads[i];  

        uint8_t hits_size;
        if (is_om) hits_size = read.hits_om.size(); 
        else hits_size = read.hits_tm.size();

        
        for (uint8_t h = 0; h < hits_size; h++){

            aligned_num_b = 0;
            aligned_num_f = 0;
            memset(aligned_sai_b, 0, sizeof(aligned_sai_b));
            memset(aligned_sai_f, 0, sizeof(aligned_sai_f));

            if (is_om){
                hit_om_t hit = read.hits_om[h]; 

                if (hit.is_aligned_bck == 1){
                    for (uint32_t j = 0; j < hit.high - hit.low && j < ALIGN_REPORT_NUM; j++){
                        aligned_sai_b[j] = sai[hit.low + j] + 1;
                        aligned_num_b++;

                        suppress_aligned += (hit.high - hit.low) > ALIGN_REPORT_NUM ? 
                            hit.high - hit.low - ALIGN_REPORT_NUM: 0;
                    }
                }
                if (hit.is_aligned_bck == 0){
                    for (uint32_t j = 0; j < hit.high - hit.low && j < ALIGN_REPORT_NUM; j++){
                        aligned_sai_f[j] = sai_rev[hit.low + j] + 1;
                        aligned_num_f++;

                        suppress_aligned += (hit.high - hit.low) > ALIGN_REPORT_NUM ? 
                            hit.high - hit.low - ALIGN_REPORT_NUM: 0;
                    }               
                }

                if ((bytes + 1024) > BUFF_SIZE) {
                    writeFile(fp_sam, buffer, bytes);
                    bytes = 0;
                }

                // write to SAM file
                for (uint32_t j = 0; j < aligned_num_b; j++){
                    writeSAMFields<hit_om_t>(bytes, 
                                            buffer,
                                            chrs,
                                            read,
                                            &hit, 
                                            suppress_aligned,
                                            aligned_sai_b,
                                            aligned_sai_f,
                                            true,
                                            j,
                                            is_om);          
                }

                for (uint32_t j = 0; j < aligned_num_f; j++){
                    writeSAMFields<hit_om_t>(bytes, 
                                            buffer,
                                            chrs,
                                            read,
                                            &hit,
                                            suppress_aligned,
                                            aligned_sai_b,
                                            aligned_sai_f,
                                            false,
                                            j,
                                            is_om);
                }
            }
            else{
                hit_tm_t hit = read.hits_tm[h]; 

                if (hit.is_rev_cmpt == 0){
                    for (uint32_t j = 0; j < hit.high - hit.low && j < ALIGN_REPORT_NUM; j++){
                        aligned_sai_b[j] = sai[hit.low + j] + 1;
                        aligned_num_b++;

                        suppress_aligned += (hit.high - hit.low) > ALIGN_REPORT_NUM ? 
                            hit.high - hit.low - ALIGN_REPORT_NUM: 0;
                    }
                }

                if (hit.is_rev_cmpt == 1){
                    for (uint32_t j = 0; j < hit.high - hit.low && j < ALIGN_REPORT_NUM; j++){
                        aligned_sai_f[j] = sai[hit.low + j] + 1;
                        aligned_num_f++;

                        suppress_aligned += (hit.high - hit.low) > ALIGN_REPORT_NUM ? 
                            hit.high - hit.low - ALIGN_REPORT_NUM: 0;
                    }               
                }

                if ((bytes + 1024) > BUFF_SIZE) {
                    writeFile(fp_sam, buffer, bytes);
                    bytes = 0;
                }

                // write to SAM file
                for (uint32_t j = 0; j < aligned_num_b; j++){
                    writeSAMFields<hit_tm_t>(bytes, 
                                            buffer,
                                            chrs,
                                            read,
                                            &hit, 
                                            suppress_aligned,
                                            aligned_sai_b,
                                            aligned_sai_f,
                                            true,
                                            j,
                                            is_om);          
                }

                for (uint32_t j = 0; j < aligned_num_f; j++){
                    writeSAMFields<hit_tm_t>(bytes, 
                                            buffer,
                                            chrs,
                                            read,
                                            &hit,
                                            suppress_aligned,
                                            aligned_sai_b,
                                            aligned_sai_f,
                                            false,
                                            j,
                                            is_om);
                }
            }
        }
    }
            
    if (bytes > 0) {
        writeFile(fp_sam, buffer, bytes);
        bytes = 0;
    }
}

#endif