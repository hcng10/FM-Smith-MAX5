//
// Created by hcng on 6/28/18.
//

#include "file_op.h"


void openFile(FILE **fp, string f_name, const char *mode) {

    *fp = fopen(f_name.c_str(), mode);
    if (!(*fp)) {
        fprintf(stderr, "error: unable to open file '%s'!\n", f_name.c_str());
        exit(1);
    }
}

uint64_t fileSizeBytes(FILE *fp){

    struct stat st;
    uint64_t len;
    int fd;

    // map a stream pointer to a file descriptor
    if ((fd = fileno(fp)) == -1) {
        printf("error: unable to get file size!\n");
        exit(1);
    }

    // get file status
    if(fstat(fd, &st)) {
        printf("error: unable to get file size!\n");
        exit(1);
    }

    // get the total size in byte
    len = st.st_size;

    return len;
}

void writeFile(FILE *fp, void *data, uint64_t n_bytes)
{
    if (fwrite(data, n_bytes, 1, fp) != 1) {
        fprintf(stderr, "error: unable to write file!");
        exit(1);
    }
}


void writeNinfo(FILE * fp, uint64_t ref_cnt, uint64_t fmt_cnt, uint32_t un_cnt, uint64_t cum_un_cnt){

    fwrite(&ref_cnt, sizeof(char), sizeof(uint64_t), fp);
    fwrite(&fmt_cnt, sizeof(char), sizeof(uint64_t), fp);
    fwrite(&un_cnt, sizeof(char), sizeof(uint32_t), fp);
    fwrite(&cum_un_cnt, sizeof(char), sizeof(uint64_t), fp);

}