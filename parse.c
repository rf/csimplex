#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DATATYPE float

int
split (char *** output, char * input, char * delim) {
  char * next;
  int count = 0;
  char ** out = NULL;
  while ((next = strsep(&input, delim)) != NULL) {
    count += 1;
    out = realloc(out, count * sizeof(char*));
    out[count - 1] = next;
  }
  *output = out;
  return count;
}

int
parse_rows (DATATYPE ** output, char * input, char * delim) {
  char ** lines;
  int i, count = split(&lines, input, delim);
  *output = malloc(count * sizeof(DATATYPE));

  for (i = 0; i < count; i++) (*output)[i] = atof(lines[i]);
  return count;
}

void
parse_matrix (int * rows, int * cols, DATATYPE *** output, char * input) {
  char ** lines;
  int i, count = split(&lines, input, "\n");
  DATATYPE ** out = malloc(count * sizeof(DATATYPE *));
  *cols = -1;

  for (i = 0; i < count; i++) {
    int num = parse_rows(&(out[i]), lines[i], " ");
    if (*cols == -1) *cols = num - 1;
  }

  *output = out;
  *rows = count - 1;
}

int
get_file (char * filename, char ** res)
{
  FILE * fp = fopen (filename, "rb");
  int size;
  if (fp == NULL)
    return -1; 
  fseek (fp, 0, SEEK_END);
  size = ftell (fp);
  fseek (fp, 0, SEEK_SET);
  *res = malloc (size + 1); 
  memset (*res, '\0', size+ 1); 
  if (size != fread (*res, sizeof (char), size, fp))
    return -2; 
  fclose (fp);
  return size;
}

