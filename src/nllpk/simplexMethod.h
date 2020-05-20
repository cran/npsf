/*
  File:             simplexMethod.h
  Created by:       Pavlo Mozharovskyi
  First published:  15.01.2018
  Last revised:     15.01.2018

  Declaration of the class "simplexMethod".

  This file is part of the "No license linear programming kit" (NlLpk).
  It can be distributed, modified and used without any restrictions.
*/

class simplexMethod{
  // Problem in the input form
  double** Ainp;
  double* Binp;
  double* Cinp;
  double* AinpRaw;
  int* Dinp;
  int AinpnRow;
  int AinpnCol;
  // Problem in the form suited for simplex method
  double** A;
  double* B;
  double* C;
  double* Araw;
  int AnRow;
  int AnCol;
  int nSlacks;
  // Tableau for Phase I of the simplex method
  double** tI;
  double* tIraw;
  int tInRow;
  int tInCol;
  int* tIbasis;
  int tInBasis;
  // Tableau for Phase I of the simplex method
  double** tII;
  double* tIIraw;
  int tIInRow;
  int tIInCol;
  int* tIIbasis;
  int tIInBasis;
  // Output
  int status;
  double optVal;
  double* coefs;
  // Scaling constants
  double abDen;
  // Settings
  int maxIter;
  double eps;
public:
  simplexMethod(){
    // Put zeros into pointers to recognize their emptiness
    A = 0;
    Araw = 0;
    B = 0;
    C = 0;
    Ainp = 0;
    AinpRaw = 0;
    Binp = 0;
    Cinp = 0;
    Dinp = 0;
    tI = 0;
    tIraw = 0;
    tIbasis = 0;
    tII = 0;
    tIIraw = 0;
    tIIbasis = 0;
    coefs = 0;
    // Initialize variables
    tIInBasis = 0;
    // Initialize constants
    maxIter = 100000;
    eps = 1e-8;
    abDen = 0;
  }
  int inputProblem(double* A, double* B, double* C, int* D,
                   int nrowA, int ncolA);
  int setProblem();
  int setTableauI();
  int printTableau(double** tableau, int nRow, int nCol);
  int printBasis(int* basis, int nBasis);
  int getPivotCol(double** tableau, int* basis, int nRow, int nCol, 
                  bool randomize);
  int getPivotRow(double** tableau, int nRow, int nCol, int pivotCol);
  int doPivot(double** tableau, int nRow, int nCol,
              int pivotCol, int pivotRow, int *basis);
  int optimize(double** tableau, int nRow, int nCol, int *basis);
  int solve();
  int setTableauII();
  int readSolution(int *status, double *optVal, double* coefs);
  int freeMemory();
  int scale();
  int backScale();
  int nullize(double** tableau, int nRow, int nCol);
};
