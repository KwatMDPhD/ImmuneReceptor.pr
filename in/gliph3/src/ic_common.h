#ifndef _COMMON_H_
#define _COMMON_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <dirent.h>
#include <signal.h>

#include "klib/kseq.h"
#include "klib/kvec.h"
#include "klib/khash.h"
#include "klib/kstring.h"
#include "klib/kmath.h"

// #define DEBUG 1
#define VERSION "3.01"
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#define AMINOACIDS "ARNDCEQGHILKMFPSTWYV"
#define CHAIN_AG "TRA"
#define CHAIN_BD "TRB"
#define NA_SYMBOL "NA"
#define NUM_HLA_ALLELE 18
#define MAXLINE 1024
#define MAX_NUM_KMER 5
#define MAX_CDR3_LENGTH 50
#define MAX_HLA_LENGTH 50
#define MAX_SID_LENGTH 100
#define H5 500
#define H2 200
#define H1 100
#define K1 1000
#define MAX_FILENAME_LENGTH 200
#define EPSILON 0.0000000001
#define NUM_OF_LINE_TO_SORT 200000
#define MOTIF_OVERLAP 2
#define ENTROPY_KMER_SIZE 3
#define ENTROPY_FRACTION_PERC 0.95
#define NA_CDR3_DISTANCE -1000
#define NA_NA_DISTANCE -500
#define OTHER_CHAIN_CUTOFF 2

#endif
