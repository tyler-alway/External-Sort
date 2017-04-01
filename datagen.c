//------------------------------------------------------------------------- 
// This is supporting software for CS415P/515 Parallel Programming.
// Copyright (c) Portland State University
//------------------------------------------------------------------------- 

// A simple data generation engine.
// 
// - With the given number N, it generates a random permutation of [1..N] 
//   and writes to stdout the byte representation of these integers
//
// Usage: 
//   linux> ./datagen N > <outfile>	-- pipe the output to a file
//
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv) 
{
  int *array;
  int N = 64;	// default value
  int i, j, tmp, intSize = sizeof(int);

  if (argc > 1)
    N = atol(argv[1]);

  array = (int *) malloc(intSize * N);
  for (i=0; i<N; i++)
    array[i] = i + 1;
  for (i=0; i<N; i++) {
    j = (rand()*1./RAND_MAX)*(N-1);
    tmp = array[i];
    array[i] = array[j];
    array[j] = tmp;
  }
  fwrite(array, intSize, N, stdout);
}

