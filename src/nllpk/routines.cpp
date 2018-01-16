/*
  File:             routines.cpp
  Created by:       Pavlo Mozharovskyi
  First published:  15.01.2018
  Last revised:     15.01.2018

  Routines for the "No license linear programming kit" (NlLpk).

  This file is part of the "No license linear programming kit" (NlLpk).
  It can be distributed, modified and used without any restrictions.
*/

#include "NlLpk.h"

// Checks whether an 'entry' is in the 'array' of length 'n'
int isInArray(int entry, int* array, int n){
  for (int i = 0; i < n; i++){
    if (entry == array[i]){
      return i;
    }
  }
  return -1;
}
