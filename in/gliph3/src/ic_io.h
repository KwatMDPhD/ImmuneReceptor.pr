#ifndef _IO_H_
#define _IO_H_

#include "ic_common.h"
#include "ic_type.h"

char* ic_readline(FILE *fp);
void ic_read_into_line_vector(const char *ifile, s_vector *array);

#endif
