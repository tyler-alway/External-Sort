//------------------------------------------------------------------------- 
// This is supporting software for CS415P/515 Parallel Programming.
// Copyright (c) Portland State University.
//------------------------------------------------------------------------- 

// An utility program for verifying that the integers in a data file
// are sorted in the ascending order.
//
// Usage: 
//   linux> ./verify <datafile>
//
//
#include <stdio.h>

int main(int argc, char **argv) 
{
  FILE *fp;
  int nbytes, val1, val2, i;
  int intSize = sizeof(int);

  if (argc != 2) {
    printf("Usage: ./verify <datafile>\n");
    return 0;
  }
  if (!(fp = fopen(argv[1], "r"))) {
    printf("Can't open datafile %s\n", argv[1]);
    return 0;
  }

  nbytes = fread(&val1, intSize, 1, fp);
  if (nbytes == 0) {
    printf("No data in %s\n", argv[1]);
    return 0;
  }
  i = 1;
  while (1) {
    nbytes = fread(&val2, intSize, 1, fp);
    if (nbytes == 0)
      break;
    i++;
    if (val1 > val2) {
      printf("Failed: items %u, %u: %u > %u\n",
	     i-1, i, val1, val2);
      return 0;
    }
    val1 = val2;
  }
  printf("Data in %s are sorted\n", argv[1]);
}

