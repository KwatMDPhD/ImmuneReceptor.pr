#include "ic_io.h"
#include "ic_util.h"

char *ic_readline(FILE *fp)
{
	int len = MAXLINE;
	char *buffer = icalloc(len, char);
	memset(buffer, 0, len);
	char c;
	int i = 0;
	do
	{
		c = fgetc(fp);
		if (feof(fp))
		{
			break;
		}
		if (i >= len - 1)
		{
			len   *= 2;
			buffer = realloc(buffer, len);
		}
		buffer[i++] = c;
		if (c == '\r' || c == '\n')
		{
			break;
		}
	} while (1);

	if (i == 0)
	{
		safe_free(buffer);
		return(NULL);
	}
	buffer[i] = '\0';
	return(buffer);
}


void ic_read_into_line_vector(const char *ifile, s_vector *array)
{
	FILE *fp = fopen(ifile, "r");

	if (fp == NULL)
	{
		fprintf(stderr, "Could not open file:%s\n", ifile);
		exit(EXIT_FAILURE);
	}

	char *line;
	while ((line = ic_readline(fp)) != NULL)
	{
		//remove \ n or \ r from end of line
		if (ic_all_space(line) == 1)
		{
			safe_free(line);
		}
		else
		{
			if(line[strlen(line) - 1] == '\r' || line[strlen(line) - 1] =='\n')
			{
				line[strlen(line) - 1] = 0;
			}
			kv_push(char *, *array, line);
		}
	}
	fclose(fp);
}
