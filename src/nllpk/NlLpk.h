/*
  File:             NlLpk.h
  Created by:       Pavlo Mozharovskyi
  First published:  15.01.2018
  Last revised:     15.01.2018

  Main header for the "No license linear programming kit" (NlLpk).

  This file is part of the "No license linear programming kit" (NlLpk).
  It can be distributed, modified and used without any restrictions.
*/

// System headers
#include <memory.h>
#include <stdlib.h>
#include <math.h>

// R-headers (used for testing purposes)
#include <Rcpp.h>
using namespace Rcpp;

// Project headers
#include "simplexMethod.h"
#include "routines.h"

// Function for doule-indexing
template<typename T> T* asMatrix(T arr, int n, int d){
  T* mat = new T[n];
  for (int i = 0; i < n; i++)
    mat[i] = arr + i*d;
  return mat;
}
