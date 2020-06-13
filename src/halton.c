/* 
 Obtain Halton Sequence
 https://stackoverflow.com/questions/42661304/implementing-4-dimensional-halton-sequence
 */

#include <R.h>
#include <Rmath.h>
#include "manyprimes.h"

double HaltonDraw(int index, int base){
 double f = 1, r = 0;
 while(index > 0){
  f = f / base;
  r = r + f* (index% base);
  index = index/base;
 }
 return r;
}

// *indices and *Hdr must be of length *N
void HaltonSeq(int *indices, int *N, int* base, double* Hdr){
 for(int i = 0; i < *N; i++){
  Hdr[i] = HaltonDraw(indices[i], *base);
 }
}

// get the primes
void Primes(int *indices, int *N, double* myprimes){
 for(int i = 0; i < *N; i++){
  myprimes[i] = my100008Primes[indices[i]];
 }
}

//
