#ifndef FORMAT_H
#define FORMAT_H

#include "file_op.h"

struct chr_t{
    string name;
    uint64_t begin;
    uint64_t end;
};


struct AlignedInfo_t{
	uint16_t seed_idx;
    uint16_t extend_idx;

    uint16_t aligned_score;

    bool operator<(const AlignedInfo_t &o) const
    {
        return aligned_score < o.aligned_score;
    }
};

struct Seed_t{
    bool reverse;
    uint8_t rdStart;

    uint8_t extendCnt;

    bool refValid[MAX_EXTD+1];
    uint64_t refStart[MAX_EXTD+1];
    uint8_t pckRef[MAX_EXTD+1][CEIL(REF_DEF_LEN * BP_BIT, 8)];

    //bool smithAligned[MAX_EXTD+1] = {false};
    //int16_t smithScore[MAX_EXTD+1] = {0};
};

struct Read_t{
    string readId;
    char seq[READ_DEF_LEN];
    char qs[READ_DEF_LEN];

    uint16_t readLen;

    uint8_t pckSeq[CEIL(READ_DEF_LEN * BP_BIT, 8)];
    uint8_t pckQS[CEIL(READ_DEF_LEN * QS_BIT, 8)];

    uint8_t pckSeq_rev[CEIL(READ_DEF_LEN * BP_BIT, 8)];
    uint8_t pckQS_rev[CEIL(READ_DEF_LEN * QS_BIT, 8)];

    bool withReverse;
    bool withNotReverse;

    Seed_t seedInfo[(READ_DEF_LEN / SEED_LEN) * 2];
    uint8_t seedCnt;
    
    uint64_t totalExtendCnt;

    priority_queue<AlignedInfo_t> q;
};



//info for kernel input
struct InToRead_t{
	uint64_t read_idx;
	uint16_t seed_idx;
    uint16_t extend_idx;
};

void faToFmt(FILE *fasta_fp, 
            char * fmt, 
            uint64_t & fmt_cnt,
            std::vector<chr_t> &chrs){

    // all the counting starts from zero
    char line[LINE_SIZE];

    bool chr_valid = false;
    chr_t tmp;
    //nchar_cluster_t nchar_cluster;

    // read FASTA lines
    while (fgets(line, LINE_SIZE, fasta_fp) != NULL) {
        // not header line
        if (line[0] != '>') {

            int i = 0;
            // parse line
            while (line[i] > 32) {

                char sym = toupper(line[i]);

                fmt[fmt_cnt] = sym;
                fmt_cnt++;
                
                i++;
            }
        }
        // get the chromosome name
        else{
            if (chr_valid){
                chrs.push_back(tmp);
            }
            tmp.name.assign(line, strlen(line));
            tmp.begin = fmt_cnt;
            chr_valid = true;
        }
        tmp.end = fmt_cnt;
    }

    if (chr_valid){
        chrs.push_back(tmp);
    }
}

inline void reset_read(Read_t & read, int value){
    memset(read.pckSeq, value, sizeof(read.pckSeq));
    memset(read.pckQS, value, sizeof(read.pckQS));

    memset(read.pckSeq_rev, value, sizeof(read.pckSeq_rev));
    memset(read.pckQS_rev, value, sizeof(read.pckQS_rev));

    read.withReverse = false;
    read.withNotReverse = false;

    memset(read.seedInfo, value, sizeof(read.seedInfo));
    read.seedCnt = 0;
    read.totalExtendCnt = 0;
    //read.smithAlignedCnt = 0;

}

inline void rdSeed(FILE *seed_fp, 
            vector<Read_t> & read_vec){

    char line[2048];
    Read_t read;
    reset_read(read, 0);
    read_vec.clear();


    int seed_cnt = 0;
    int extend_cnt = 0;

    // get read ID
    char * fgets_rd = fgets(line, 2048, seed_fp);
    if (fgets_rd != NULL) line[strlen(line)-1] = '\0';
    else return;

    //FILE * write_fp = NULL;

    for (int i=0; i < MAX_RD_CNT; i++){

        seed_cnt = 0;

        if (fgets_rd == NULL){
            break;
        }
        
        read.readId.assign(line, strlen(line));

        // get the sequence
        fgets_rd = fgets(line, 2048, seed_fp);
        line[strlen(line)-1] = '\0';
        strcpy(read.seq, line);
        read.readLen = strlen(line);

        // get the +
        fgets_rd = fgets(line, 2048, seed_fp);

        // get QS
        fgets_rd = fgets(line, 2048, seed_fp);
        line[strlen(line)-1] = '\0';
        strcpy(read.qs, line);

        // get all the seeds info for extension
        while(true){
            
            // get next line of seed for extension (with position so that I can return to)
            unsigned long position;
            fflush(seed_fp);
            position = ftell(seed_fp);

            fgets_rd = fgets(line, 2048, seed_fp);
            // Mark the null char
            if (fgets_rd != NULL ) line[strlen(line)-1] = '\0';            
            
            // it turns out to be next read
            if (fgets_rd == NULL ||  line[0] == '@'){
                // a new object will be created before pushing to vector
                read_vec.push_back(read);
                
                
                reset_read(read, 0);

                strcpy(read.seq, line);

                if ( i == MAX_RD_CNT-1){
                    fseek(seed_fp, position, SEEK_SET);
                }

                break;
            }
    
            Seed_t & seed = read.seedInfo[seed_cnt];

            extend_cnt = 0;

            // get the start at read
            char * pch = strtok(line, " ");
            seed.rdStart = (uint8_t) atoi(pch);

            // get if it is reversed
            pch = strtok (NULL, " ");
            if (pch[0] == 'f'){
                seed.reverse = true;
                read.withReverse = true;
            }else{
                seed.reverse = false;
                read.withNotReverse= true;
            }
            pch = strtok (NULL, " ");

            seed.extendCnt = 0;

                // get all the locations in reference
            while (pch != NULL){

                if (extend_cnt > MAX_EXTD){
                    cerr<<"More than "<<MAX_EXTD<<" extension\n";
                }
                seed.refStart[extend_cnt] = (uint64_t) atol(pch);
                seed.refValid[extend_cnt] = true;
                
                extend_cnt++;
                seed.extendCnt++;

                pch = strtok (NULL, " ");
            }

            seed_cnt++;
            read.seedCnt++;               
        }
    }

    return;
}



inline void pack_bp(uint8_t *bp, uint8_t val, uint8_t i){

	uint32_t bp_i = (i * BP_BIT)/ 8;
	// within the 8-bit space or outside the 8-bit space
	uint32_t cur_mod = (i * BP_BIT)% 8;
	uint32_t next_mod = ((i + 1) * BP_BIT)% 8;

    uint8_t val_sft;
    uint8_t next_val_sft;

    
	if (cur_mod < next_mod || next_mod == 0){
		val_sft = val << cur_mod;
		bp[bp_i] |=  val_sft;
	}else{
		val_sft = val << cur_mod;
		next_val_sft = val >> (8 - cur_mod);

		bp[bp_i] |= val_sft;
		bp[bp_i+1] |=  next_val_sft;
	}

}

inline void pck_qs(uint8_t *pckqs, char cur_qs, uint8_t i){

	uint32_t qs_i = (i * QS_BIT)/ 8;
	// within the 8-bit space or outside the 8-bit space
	uint32_t cur_mod = (i * QS_BIT)% 8;
	uint32_t next_mod = ((i + 1) * QS_BIT)% 8;

    uint8_t qs_sft;
    uint8_t next_qs_sft;

    cur_qs = cur_qs - 33;

    //cout<<"mod "<< (int)cur_mod<<" "<<(int)next_mod<<"\n\n";

	if (cur_mod < next_mod || next_mod == 0){
		qs_sft = cur_qs << cur_mod;
		pckqs[qs_i] |=  qs_sft;
	}else{
		qs_sft = cur_qs << cur_mod;
		next_qs_sft = cur_qs >> (8 - cur_mod);

		pckqs[qs_i] |= qs_sft;
		pckqs[qs_i+1] |=  next_qs_sft;
	}
}

// pack the sequence and QS
inline void pck_seq_qs(char *seq, uint8_t * pck, 
                    char *qs, uint8_t * pckqs,
                    uint8_t len){

    for (uint8_t i = 0; i < len; i++) {
    	switch(seq[i]){	
            case 'A': pack_bp(pck, 0, i); break;
            case 'C': pack_bp(pck, 1, i); break;
            case 'G': pack_bp(pck, 2, i); break;
            case 'T': pack_bp(pck, 3, i); break;
            case 'N': pack_bp(pck, 4, i); break;
            case 'a': pack_bp(pck, 0, i); break;
            case 'c': pack_bp(pck, 1, i); break;
            case 'g': pack_bp(pck, 2, i); break;
            case 't': pack_bp(pck, 3, i); break;
            case 'n': pack_bp(pck, 4, i); break;
            default: pack_bp(pck, 5, i); break;
		}

        pck_qs(pckqs, qs[i], i);
    }
}



//pack the reference
inline void pck_ref(Read_t & read, Seed_t & seed, char * ref_fmt, uint64_t ref_fmt_len){

    uint8_t rdStart = seed.rdStart;

    for (uint8_t s = 0; s < seed.extendCnt; s++){

        uint64_t refStart = seed.refStart[s];
        
        //cout<<"ho dai jei "<<refStart<<"\n";

        if (refStart - rdStart + REF_DEF_LEN >= ref_fmt_len){
            seed.refValid[s] = false;
        }else{
         
            for (uint64_t i = 0; i < REF_DEF_LEN; i++) {
                switch(ref_fmt[i + refStart - rdStart]){	
                    case 'A': pack_bp(seed.pckRef[s], 0, i); break;
                    case 'C': pack_bp(seed.pckRef[s], 1, i); break;
                    case 'G': pack_bp(seed.pckRef[s], 2, i); break;
                    case 'T': pack_bp(seed.pckRef[s], 3, i); break;
                    case 'N': pack_bp(seed.pckRef[s], 4, i); break;
                    default: pack_bp(seed.pckRef[s], 5, i); break;            
                }
            }//cout<<"\n";
            read.totalExtendCnt = read.totalExtendCnt + 1;
        }

    }    
}




inline void packRdRef(vector<Read_t> & read_vec, 
                vector<InToRead_t> & intoread_vec, 
                char * ref_fmt, 
                uint64_t ref_len, 
                uint64_t & total_extend_cnt){
#if IS_SIM == 1
    cout<<"The vector size is....."<<read_vec.size()<<" \n";
#endif

    intoread_vec.clear();
    if (read_vec.size() == 0) return;

    int n_threads = omp_get_max_threads();    
#pragma omp parallel for num_threads(n_threads)
    for (uint i = 0; i < read_vec.size(); i++){
    
        Read_t &cur_read = read_vec.at(i);  
        uint8_t read_len = strlen(cur_read.seq);

        //1. pack read, and qs for the normal one
        if (cur_read.withNotReverse){
            pck_seq_qs(cur_read.seq, cur_read.pckSeq, 
                    cur_read.qs, cur_read.pckQS, 
                    cur_read.readLen);
        }
        
        //2. pack read, and qs for the reverse complement
        if (cur_read.withReverse){
            char seq_rev[READ_DEF_LEN];
            char qs_rev[READ_DEF_LEN];

            strcpy(seq_rev, cur_read.seq);
            strcpy(qs_rev, cur_read.qs);

            for (size_t s = 0; s < (size_t) (cur_read.readLen / 2); s++){
                swap(seq_rev[s], seq_rev[cur_read.readLen-1-s]);
                swap(qs_rev[s], qs_rev[cur_read.readLen-1-s]);
            }

            for (size_t s = 0; s < cur_read.readLen; s++){
                // after reverse, we make complement
                char c = 'N';
                switch(seq_rev[s]){
                    case 'A': c = 'T'; break;
                    case 'T': c = 'A'; break;
                    case 'C': c = 'G'; break;
                    case 'G': c = 'C'; break;
                    case 'a': c = 'T'; break;
                    case 't': c = 'A'; break;
                    case 'c': c = 'G'; break;
                    case 'g': c = 'C'; break;
                    case 'N': c = 'N'; break;
                }
                seq_rev[s] = c;
            }

#if IS_SIM == 1
            cout<<"original: "<<cur_read.seq<<"\n";
            cout<<"reverse: "<<seq_rev<<"\n";

            cout<<"original: "<<cur_read.qs<<"\n";
            cout<<"reverse: "<<qs_rev<<"\n";
#endif

            pck_seq_qs(seq_rev, cur_read.pckSeq_rev, 
                    qs_rev, cur_read.pckQS_rev, 
                    cur_read.readLen);

            //for (int gg =0;gg<CEIL(READ_DEF_LEN * BP_BIT, 8);gg++){
                //std::cout << gg<<"seq = " << std::bitset<8>(cur_read.pckSeq_rev[gg])  << std::endl;
            //}
        }

        for (uint s = 0; s < cur_read.seedCnt; s++){
            pck_ref(cur_read, cur_read.seedInfo[s], ref_fmt, ref_len);
        }        
    }

    //counting the total number of extends in the end
    //record the extend info
    for (uint64_t i = 0; i < read_vec.size(); i++){
        total_extend_cnt = total_extend_cnt + read_vec[i].totalExtendCnt; 

        for (uint s = 0; s < read_vec[i].seedCnt; s++){
            for (uint e = 0; e < read_vec[i].seedInfo[s].extendCnt; e++){

                InToRead_t tmp_info;

                tmp_info.read_idx = i;
                tmp_info.seed_idx = (uint16_t) s;
                tmp_info.extend_idx = (uint16_t) e;

                if (read_vec[i].seedInfo[s].refValid[e]  == true){
                    intoread_vec.push_back(tmp_info);

                    #if IS_SIM == 1
                    cout<< i<<" "<<s <<" "<<e<<" valid "<<read_vec[i].seedInfo[s].refValid[e]<<" "<<(int)read_vec[i].seedInfo[s].extendCnt<<"\n";
                    #endif
                }

                
            }
        }    
    }

    return;
}


/*
inline void rdSeed(FILE *seed_fp, 
            vector<Read_t> & read_vec){

    char line[2048];
    Read_t read;
    reset_read(read, 0);


    int seed_cnt = 0;
    int extend_cnt = 0;

    // get read ID
    char * fgets_rd = fgets(line, 2048, seed_fp);
    if (fgets_rd != NULL) line[strlen(line)-1] = '\0';

    //FILE * write_fp = NULL;
    for (int iter_i = 0; ;iter_i++){
        if (fgets_rd == NULL){
            break;
        }

        for (int i=0; i < MAX_RD_CNT; i++){

            seed_cnt = 0;
            read.withReverse = false;
            read.withNotReverse = false;
            
            read.seedCnt = 0;
            read.smithAlignedCnt = 0;


            if (fgets_rd == NULL){
                break;
            }
            read.readId.assign(line, strlen(line));

            // get the sequence
            fgets_rd = fgets(line, 2048, seed_fp);
            line[strlen(line)-1] = '\0';
            strcpy(read.seq, line);
            read.readLen = strlen(line);

            // get the +
            fgets_rd = fgets(line, 2048, seed_fp);

            // get QS
            fgets_rd = fgets(line, 2048, seed_fp);
            line[strlen(line)-1] = '\0';
            strcpy(read.qs, line);

            // get all the seeds info for extension
            while(true){
                
                // get next line of seed for extension
                fgets_rd = fgets(line, 2048, seed_fp);
                // Mark the null char
                if (fgets_rd != NULL ) line[strlen(line)-1] = '\0';            
                
                // it turns out to be next read
                if (fgets_rd == NULL ||  line[0] == '@'){
                    // a new object will be created before pushing to vector
                    read_vec.push_back(read);

                    
                    reset_read(read, 0);

                    strcpy(read.seq, line);

                    break;
                }
      
                Seed_t & seed = read.seedInfo[seed_cnt];

                extend_cnt = 0;

                // get the start at read
                char * pch = strtok(line, " ");
                seed.rdStart = (uint8_t) atoi(pch);

                // get if it is reversed
                pch = strtok (NULL, " ");
                if (pch[0] == 'f'){
                    seed.reverse = true;
                    read.withReverse = true;
                }else{
                    seed.reverse = false;
                    read.withNotReverse= true;
                }
                pch = strtok (NULL, " ");

                seed.extendCnt = 0;

                 // get all the locations in reference
                while (pch != NULL){

                    if (extend_cnt > MAX_EXTD){
                        cerr<<"More than "<<MAX_EXTD<<" extension\n";
                    }
                    seed.refStart[extend_cnt] = (uint64_t) atol(pch);
                    seed.refValid[extend_cnt] = true;
                    
                    extend_cnt++;
                    seed.extendCnt++;

                    pch = strtok (NULL, " ");
                }

                seed_cnt++;
                read.seedCnt++;               
            }
        }
    }

}
*/

#endif 
