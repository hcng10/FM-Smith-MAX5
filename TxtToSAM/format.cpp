#include "format.h"


inline void setSLen(uint64_t & s_position, uint64_t & s_len, uint64_t & i){
    s_position = i;
	s_len = 0;
}

void loadtxt(FILE *fp, 
                std::vector<read2Bit_t> &reads2Bit, 
                bool is_om,
                char *buffer, 
                uint64_t size_r, 
                uint64_t *bytes, 
                bool r_ctrl){
        
    reads2Bit.clear();

    if (r_ctrl == true && size_r > 0) {
		// read a bunch of things from file first
		if (fread(buffer, size_r, 1, fp) != 1) {
			fprintf(stderr, "error: unable to read file!\n");
			exit(1);
		}

        // Then read until new line
        int c;
        while ((c = fgetc(fp)) != EOF) {
            buffer[size_r++] = (char)c;
            if (c == '\n') {
                break;
            }
            //buffer[size_r++] = (char)c;
        }

        // update bytes read
        *bytes += size_r;

        char char_buf[100];
        //string tmp_at_line;

        // parse the buffer
        read2Bit_t tmp2Bit;
        hit_om_t tmpHit_om;
        hit_tm_t tmpHit_tm;

        uint64_t s_position;
        uint64_t s_len;
        uint64_t i = 0;

        
        while (i < size_r) {
            tmp2Bit.hits_om.clear();
            tmp2Bit.hits_tm.clear();

            //// Step 1: get meta data line, which is <@r337>
            setSLen(s_position, s_len, i);

            while (buffer[i++] != ' ') {
                s_len += 1;
            }
            tmp2Bit.at_line.assign(&buffer[s_position], s_len);

            //// Step 2: get seq len
            setSLen(s_position, s_len, i);

            while (buffer[i++] != ' ') {
                s_len += 1;
            }
            tmp2Bit.seq.assign(&buffer[s_position], s_len);
            tmp2Bit.seq_len = s_len;

            //// Step 3: get q_score
            setSLen(s_position, s_len, i);

            while (buffer[i++] != ' ') {
                s_len += 1;
            }
            tmp2Bit.q_score.assign(&buffer[s_position], s_len);
        
            ////Step 4: get the number of hits
            //bool unknown = false;
            setSLen(s_position, s_len, i);

            while (buffer[i++] != ' ') {
                s_len += 1;
            }
            memcpy(char_buf, &buffer[s_position], s_len);
            char_buf[s_len] = '\0';
            tmp2Bit.n_hits = atoi(char_buf);

            memset(char_buf, 0, s_len+1);
            setSLen(s_position, s_len, i);
    //cout<<tmp2Bit.at_line<<" "<<(int)tmp2Bit.seq_len<<"\n";
            for (uint8_t h = 0; h < tmp2Bit.n_hits; h++){

                if (is_om){           
                    // rev_cmpt?
                    tmpHit_om.is_rev_cmpt = (buffer[i++] == '1') ? true: false;
                    i++;
                    setSLen(s_position, s_len, i);

                    // backword?
                    tmpHit_om.is_aligned_bck = (buffer[i++] == '1') ? true: false;
                    i++;
                    setSLen(s_position, s_len, i);

                    // low
                    while (buffer[i++] != ' ') {
                        s_len += 1;
                    }
                    memcpy(char_buf, &buffer[s_position], s_len);
                    char_buf[s_len] = '\0';
                    tmpHit_om.low = strtoul(char_buf, NULL, 0);

                    memset(char_buf, 0, s_len+1);
                    setSLen(s_position, s_len, i);

                    //high
                    while (buffer[i++] != ' ') {
                        s_len += 1;
                    }
                    memcpy(char_buf, &buffer[s_position], s_len);
                    char_buf[s_len] = '\0';
                    tmpHit_om.high = strtoul(char_buf, NULL, 0);

                    memset(char_buf, 0, s_len+1);
                    setSLen(s_position, s_len, i); 

                    //mispos
                    while (buffer[i++] != ' ') {
                        s_len += 1;
                    }
                    memcpy(char_buf, &buffer[s_position], s_len);
                    char_buf[s_len] = '\0';
                    tmpHit_om.mis_pos = (uint8_t) strtoul(char_buf, NULL, 0);

                    memset(char_buf, 0, s_len+1);
                    setSLen(s_position, s_len, i); 

                    // missym
                    tmpHit_om.mis_sym = buffer[i++];
                    i++;
                    setSLen(s_position, s_len, i);

                    // qs
                    tmpHit_om.qs = buffer[i++];
                    i++;
                    setSLen(s_position, s_len, i);

                    tmp2Bit.hits_om.push_back(tmpHit_om);
                    //cout<<s_position <<" "<<s_len<<" "<<i<<"\n";
                    //cout<< tmpHit_om.is_rev_cmpt <<" "<<tmpHit_om.low<<" "<<tmpHit_om.high <<" "<<(int)tmpHit_om.mis_pos <<" "<<tmpHit_om.mis_sym<<" "<<tmpHit_om.qs<<" "<<i<<" "<<size_r<<"\n";
                }
                else{
                    // rev_cmpt?
                    tmpHit_tm.is_rev_cmpt = (buffer[i++] == '1') ? true: false;
                    i++;
                    setSLen(s_position, s_len, i);

                    // low
                    while (buffer[i++] != ' ') {
                        s_len += 1;
                    }
                    memcpy(char_buf, &buffer[s_position], s_len);
                    char_buf[s_len] = '\0';
                    tmpHit_tm.low = strtoul(char_buf, NULL, 0);

                    memset(char_buf, 0, s_len+1);
                    setSLen(s_position, s_len, i);

                    //high
                    while (buffer[i++] != ' ') {
                        s_len += 1;
                    }
                    memcpy(char_buf, &buffer[s_position], s_len);
                    char_buf[s_len] = '\0';
                    tmpHit_tm.high = strtoul(char_buf, NULL, 0);

                    memset(char_buf, 0, s_len+1);
                    setSLen(s_position, s_len, i); 

                    //mispos0
                    while (buffer[i++] != ' ') {
                        s_len += 1;
                    }
                    memcpy(char_buf, &buffer[s_position], s_len);
                    char_buf[s_len] = '\0';
                    tmpHit_tm.mis_pos0 = (uint8_t) strtoul(char_buf, NULL, 0);

                    memset(char_buf, 0, s_len+1);
                    setSLen(s_position, s_len, i); 

                    // missym0
                    tmpHit_tm.mis_sym0 = buffer[i++]; i++;
                    // qs0
                    tmpHit_tm.qs0 = buffer[i++]; i++;
                    setSLen(s_position, s_len, i);


                    //mispos1
                    while (buffer[i++] != ' ') {
                        s_len += 1;
                    }
                    memcpy(char_buf, &buffer[s_position], s_len);
                    char_buf[s_len] = '\0';
                    tmpHit_tm.mis_pos1 = (uint8_t) strtoul(char_buf, NULL, 0);

                    memset(char_buf, 0, s_len+1);
                    setSLen(s_position, s_len, i); 

                    // missym1
                    tmpHit_tm.mis_sym1 = buffer[i++]; i++;
                    // qs1
                    tmpHit_tm.qs1 = buffer[i++]; i++;
                    setSLen(s_position, s_len, i);

                    tmp2Bit.hits_tm.push_back(tmpHit_tm);

                    //cout<< tmpHit_tm.is_rev_cmpt <<" "<<tmpHit_tm.low<<" "<<tmpHit_tm.high <<" "
                    //<< (int)tmpHit_tm.mis_pos0 <<" "<<tmpHit_tm.mis_sym0<<" "<<tmpHit_tm.qs0<<" "
                    //<< (int)tmpHit_tm.mis_pos1 <<" "<<tmpHit_tm.mis_sym1<<" "<<tmpHit_tm.qs1<<"\n";//" "<<i<<" "<<size_r
                }
            }
            i++;
            setSLen(s_position, s_len, i);//cout<<"\n";

            reads2Bit.push_back(tmp2Bit);
        } 
    }
}

void readChrName(FILE * fp, uint16_t chrs_num, std::vector<chr_t> &chrs){

    chr_t chr;
    uint64_t length;
    char name[NAME_LEN];

    for (uint16_t i = 0; i < chrs_num; i++){
        fread (&(chr.begin), sizeof(uint64_t), 1, fp);
        fread(&(chr.end), sizeof(uint64_t), 1, fp);
        fread(&length, sizeof(uint64_t), 1, fp);
        fread(&name, sizeof(char), length, fp);

        name[length] = '\0';

        //chr.name = name;
        chr.name = string(name, length); 

        chrs.push_back(chr);
    }
}