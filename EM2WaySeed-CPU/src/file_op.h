/*
 * file_op.h
 *
 *  Created on: Jan 31, 2019
 *      Author: hn915
 */

#ifndef FILE_OP_H_
#define FILE_OP_H_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>

#include "fm-sa-def.h"

using namespace std;

// open file
void openFile(FILE **fp, string f_name, const char *mode);
void readFile(FILE *fp, void *a, uint64_t n_bytes);
uint64_t fileSizeBytes(FILE *fp);
void writeFile(FILE *fp, void *a, uint64_t n_bytes);


#endif /* FILE_OP_H_ */
