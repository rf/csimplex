#pragma once

void parse_matrix (int * rows, int * cols, DATATYPE *** output, char * input);
int parse_rows (DATATYPE ** output, char * input, char * delim);
int split (char *** output, char * input, char * delim);
int get_file (char * filename, char ** res);
