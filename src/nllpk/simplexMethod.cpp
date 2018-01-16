/*
  File:             simplexMethod.cpp
  Created by:       Pavlo Mozharovskyi
  First published:  15.01.2018
  Last revised:     15.01.2018
 
  Implementation of the class "simplexMethod".
 
  This file is part of the "No license linear programming kit" (NlLpk).
  It can be distributed, modified and used without any restrictions.
*/

#include "NlLpk.h"

// Copy the input problem into the preliminary structure
int simplexMethod::inputProblem(double* A, double* B, double* C, int* D,
                 int nrowA, int ncolA){
  // Clean the memory if occupied
  if (this->Ainp){
    delete[] this->Ainp;
    delete[] this->AinpRaw;
  }
  if (this->Binp){
    delete[] this->Binp;
  }
  if (this->Cinp){
    delete[] this->Cinp;
  }
  if (this->Dinp){
    delete[] this->Dinp;
  }
  // Determine size and allocate structures
  this->AinpnRow = nrowA;
  this->AinpnCol = ncolA;
  this->AinpRaw = new double[this->AinpnRow * this->AinpnCol];
  this->Ainp = asMatrix(AinpRaw, this->AinpnRow, this->AinpnCol); // new double*
  this->Binp = new double[this->AinpnRow];
  this->Cinp = new double[this->AinpnCol];
  this->Dinp = new int[this->AinpnRow];
  // Copy the structures
  memcpy(this->AinpRaw, A, this->AinpnRow * this->AinpnCol * sizeof(double));
  memcpy(this->Binp, B, this->AinpnRow * sizeof(double));
  memcpy(this->Cinp, C, this->AinpnCol * sizeof(double));
  memcpy(this->Dinp, D, this->AinpnRow * sizeof(int));
//  Rcout << "Constraints: ";
//  for (int i = 0; i < AinpnRow; i++){
//    Rcout << Dinp[i] << " ";
//  }
//  Rcout << std::endl;
  return 0;
}

// Prepare program in the simplex method structure
int simplexMethod::setProblem(){
  // Clean the memory if occupied
  if (A){
    delete[] A;
    delete[] Araw;
  }
  if (B){
    delete[] B;
  }
  if (C){
    delete[] C;
  }
  // Determine number of slack variables to be added
  nSlacks = 0;
  for (int i = 0; i < AinpnRow; i++){
    if (Dinp[i] > 0){
      nSlacks++;
    }
  }
//  Rcout << nSlacks << " slacks introduced.\n";
  // Determine size and allocate structures
  AnRow = AinpnRow;
  AnCol = AinpnCol + nSlacks;
  Araw = new double[AnRow * AnCol];
  A = asMatrix(Araw, AnRow, AnCol); // new double*
  B = new double[AnRow];
  C = new double[AnCol];
//  Rcout << "AnRow=" << AnRow << ",AnCol=" << AnCol << std::endl;
  // Copy and augment the structures:
  // Constraints
  int iSlack = 0;
  for (int i = 0; i < AnRow; i++){
    if (Dinp[i] == 0){ // if "="
      // Just copy and fill slacks with zeros
      memcpy(A[i], Ainp[i], AinpnCol * sizeof(double));
      for (int j = AinpnCol; j < AnCol; j++){
        A[i][j] = 0; // no slack variables
      }
      B[i] = Binp[i]; // copy right-hand side
    }
    if (Dinp[i] == 1){ // if "<="
      // Copy and add slack
      memcpy(A[i], Ainp[i], AinpnCol * sizeof(double));
      for (int j = AinpnCol; j < AnCol; j++){
        if (j - AinpnCol == iSlack){
          A[i][j] = 1; // add slack variable where needed
        }else{
          A[i][j] = 0;
        }
      }
      iSlack++;
      B[i] = Binp[i]; // copy right-hand side
    }
    if (Dinp[i] == 2){ // if ">="
      // Copy multiplied with -1 and add slack
      for (int j = 0; j < AinpnCol; j++){
        A[i][j] = -Ainp[i][j];
      }
      for (int j = AinpnCol; j < AnCol; j++){
        if (j - AinpnCol == iSlack){
          A[i][j] = 1; // add slack variable where needed
        }else{
          A[i][j] = 0;
        }
      }
      iSlack++;
      B[i] = -Binp[i]; // copy negated right-hand side
    }
  }
  // Objective function
  memcpy(C, Cinp, AinpnCol * sizeof(double));
  for (int i = AinpnCol; i < AnCol; i++){
    C[i] = 0; // for slack variables
  }
  return 0;
}

// Create the tableau for Phase I of the simplex method
int simplexMethod::setTableauI(){
  // Clean the memory if occupied
  if (tI){
    delete[] tI;
    delete[] tIraw;
  }
  // Determine size and allocate structures
  tInRow = AnRow;
  tInCol = AnCol + AnRow;
  tIraw = new double[(tInRow + 1) * (tInCol + 1)];
  tI = asMatrix(tIraw, tInRow + 1, tInCol + 1); // new
  // Set zero row to zeros
  for (int i = 0; i <= tInCol; i++){
    tI[0][i] = 0;
  }
  // Fill artificial variables part with zeros/ones
  for (int i = 1; i <= tInRow; i++){
    for (int j = tInRow + nSlacks + 1; j <= tInCol; j++){
      if (j == i + tInRow + nSlacks){
        tI[i][j] = 1;
      }else{
        tI[i][j] = 0;
      }
    }
  }
  // Fill the basis structure
  tInBasis = tInRow;
  tIbasis = new int[tInBasis];
  for (int i = 0; i < tInRow; i++){
    tIbasis[i] = tInRow + nSlacks + 1 + i;
  }
  // Fill the tableau
  for (int i = 0; i < AnRow; i++){
    if (B[i] > 0){ // For a positive constraint
      // Fill the B-entry
      tI[i + 1][0] = B[i];
      tI[0][0] -= B[i];
      // Fill the A-entries
      for (int j = 0; j < AnCol; j++){
        tI[i + 1][j + 1] = A[i][j];
        tI[0][j + 1] -= A[i][j];
      }
    }else{ // For a negative constraint
      // Fill the B-entry
      tI[i + 1][0] = -B[i];
      tI[0][0] += B[i];
      // Fill the A-entries
      for (int j = 0; j < AnCol; j++){
        tI[i + 1][j + 1] = -A[i][j];
        tI[0][j + 1] += A[i][j];
      }
    }
  }
  return 0;
}

// Print the simplex tableau
int simplexMethod::printTableau(double** tableau, int nRow, int nCol){
  for (int i = 1; i <= nCol; i++){
//    Rcout << round(tableau[0][i] * 10000) / 10000 << "\t";
  }
//  Rcout << "|" << round(tableau[0][0] * 10000) / 10000 << std::endl;
  for (int i = 1; i <= nCol; i++){
//    Rcout << "--------";
  }
//  Rcout << "---------" << std::endl;
  for (int i = 1; i <= nRow; i++){
    for (int j = 1; j <= nCol; j++){
//      Rcout << round(tableau[i][j] * 10000) / 10000 << "\t";
    }
//    Rcout << "|" << round(tableau[i][0] * 10000) / 10000 << std::endl;
  }
//  Rcout << std::endl;
  return 0;
}

int simplexMethod::printBasis(int* basis, int nBasis){
//  Rcout << "Current basis: \t";
  for (int i = 0; i < nBasis; i++){
//    Rcout << basis[i] << "\t";
  }
//  Rcout << std::endl;
//  Rcout << std::endl;
  return 0;
}

// Pick a non-basis column having negative goal coefficient
int simplexMethod::getPivotCol(double** tableau, int* basis, int nRow, int nCol, 
                               bool randomize){
  // Collect non-basis columns with negative coefficients
  int* negNonBasis = new int[nCol - nRow];
  int nNegNonBasis = 0;
  for (int i = 1; i <= nCol; i++){
    // If not in basis and negative
    if (isInArray(i, basis, nRow) < 0 && tableau[0][i] < -eps){
      // Consider for picking
//      Rcout << "Pick " << i << " ";
      negNonBasis[nNegNonBasis++] = i;
    }
//    Rcout << std::endl;
  }
  // Check whether the solution is already found
  if (nNegNonBasis == 0){
    delete[] negNonBasis;
    return 0;
  }
  // Otherwise, pick one column randomly
  int iPivotCol = -1;
  if (nNegNonBasis == 1){
    iPivotCol = negNonBasis[0];
  }else{
    //int ri = rand() % nNegNonBasis; // use randomization to break inf loop
    if (randomize){
      int ri = unif_rand() * nNegNonBasis; // use randomization to break inf loop
      iPivotCol = negNonBasis[ri];
//      Rcout << "Randomizing" << std::endl;
    }else{
      iPivotCol = negNonBasis[nNegNonBasis - 1];
    }
  }
  delete[] negNonBasis;
  return iPivotCol;
}

int simplexMethod::getPivotRow(double** tableau, int nRow, int nCol, int pivotCol){
  // Create the temporary structures
  int iPivot = -1; // pivot row
  double minRatio = DBL_MAX; // minimum ratio
  bool posFound = false; // at least one entry in the column is positive
  for (int i = 1; i <= nRow; i++){
    // Check boundness
    if (tableau[i][pivotCol] <= eps){
      continue;
    }else{
      posFound = true;
      // Compute the ratio
      double tmpRatio = tableau[i][0] / tableau[i][pivotCol];
      // Compare it to the smallest until now and update if necessary
      if (tmpRatio < minRatio){
        iPivot = i;
        minRatio = tmpRatio;
      }
    }
  }
  if (posFound){
    // Return the index of the pivot row
    return iPivot;
  }else{
    // Problem unbounded
    return -1;
  }
}

int simplexMethod::doPivot(double** tableau, int nRow, int nCol,
                           int pivotCol, int pivotRow, int *basis){
  // For the pivot row, divide through the pivot element
  double mult = tableau[pivotRow][pivotCol];
  for (int j = 0; j <= nCol; j++){
    tableau[pivotRow][j] /= mult;
  }
  // Go through all rows and columns
  for (int i = 0; i <= nRow; i++){
    // If pivot row - skip, it has been just divided before
    if (i == pivotRow){
      continue;
    }
    // Determine row multiplier
    double mult = tableau[i][pivotCol] / tableau[pivotRow][pivotCol];
    // Do row operation, i.e. subtract pivot row multiplied with 'mult'
    for (int j = 0; j <= nCol; j++){
      tableau[i][j] -= tableau[pivotRow][j] * mult;
    }
  }
  // Change the basis record
  basis[pivotRow - 1] = pivotCol;
  return 0;
}

int simplexMethod::nullize(double** tableau, int nRow, int nCol){
  for (int i = 0; i <= nRow; i++){
    for (int j = 0; j <= nCol; j++){
      if (fabs(tableau[i][j]) <= pow(eps, 2)){
        tableau[i][j] = 0;
      }
    }
  }
  return 0;
}

int simplexMethod::optimize(double** tableau, int nRow, int nCol, int *basis){
  int iter = 0;
  while (iter < maxIter){
    iter++;
    int iPivotCol = -1;
    if (iter % (nRow > nCol ? nRow : nCol) == 0){
      // If we are looping for quite some time - choose a pivot column at random
      iPivotCol = getPivotCol(tableau, basis, nRow, nCol, true);
    }else{
      iPivotCol = getPivotCol(tableau, basis, nRow, nCol, false);
    }
    if (iPivotCol == 0){
//      Rcout << "Optimal solution found in " << iter << " iterations."<< std::endl;
//      printTableau(tableau, nRow, nCol);
//      printBasis(basis, nRow);
      return 0; // optimal solution found
    }else{
//      Rcout << "Pivot col: " << iPivotCol << std::endl;
      int iPivotRow = getPivotRow(tableau, nRow, nCol, iPivotCol);
//      Rcout << "Pivot row: " << iPivotRow << std::endl;
//      Rcout << std::endl;
      if (iPivotRow < 0){
//        Rcout << "Problem unbounded!" << std::endl;
//        printTableau(tableau, nRow, nCol);
//        printBasis(basis, nRow);
        return 1; // problem unbounded
      }else{
        doPivot(tableau, nRow, nCol, iPivotCol, iPivotRow, basis);
        nullize(tableau, nRow, nCol);
      }
    }
//    printTableau(tableau, nRow, nCol);
//    printBasis(basis, nRow);
  }
//  Rcout << "Maximum number of iterations exceeded!" << std::endl;
//  printTableau(tableau, nRow, nCol);
//  printBasis(basis, nRow);
  return 2; // maximum number of iterations exceeded
}

int simplexMethod::solve(){
  GetRNGstate();
  scale();
  setProblem();
  setTableauI();
//  nullize(tI, tInRow, tInCol);
//  printTableau(tI, tInRow, tInCol);
//  printBasis(tIbasis, tInRow);
  status = optimize(tI, tInRow, tInCol, tIbasis);
  if (status != 0 || fabs(tI[0][0]) > eps){
//    Rcout << "Initial feasible solution not found" << std::endl;
    if (status == 0){
      status = 3;
    }
    PutRNGstate();
    return 0;
  }
  setTableauII();
//  nullize(tII, tIInRow, tIInCol);
//  printTableau(tII, tIInRow, tIInCol);
//  printBasis(tIIbasis, tIInRow);
  status += optimize(tII, tIInRow, tIInCol, tIIbasis) * 10;
  backScale();
  // Fill the solution
  optVal = -tII[0][0]; // read the optimal value
  // Create structure for coefficients
  if (coefs){
    delete[] coefs;
  }
  coefs = new double[AinpnCol];
  // Fill structure for coefficients with zeros
  for (int i = 0; i < AinpnCol; i++){
    coefs[i] = 0;
  }
  // Fill non-zero coefficients
  for (int i = 0; i < tIInBasis; i++){
    if (tIIbasis[i] <= AinpnCol){
      coefs[tIIbasis[i] - 1] = tII[i + 1][0];
    }
  }
  PutRNGstate();
  return 0;
}

// Construct the tableau for Phase II from this for Phase I by driving
// artificial variables out of the basis
int simplexMethod::setTableauII(){
  // 1) Drive artificial variables with at least one non-zero basis entry
  // out of the basis and count those with all zero basis entries
  int nZeroBasis = 0;
  for (int i = 0; i < tInBasis; i++){
    if (tIbasis[i] > AnCol){ // if artificial variable
//      Rcout << "Driving var " << tIbasis[i] << " out of basis." << std::endl;
      // Search for a non-zero basis entry
      int iPivotCol = 0;
      for(int j = 1; j <= AnCol; j++){ // basis entries of this row
        if (fabs(tI[i + 1][j]) > eps){ // if non-zero entry
          iPivotCol = j; // choose this column as pivot
          break;
        }
      }
//      Rcout << "Pivot col: " << iPivotCol << std::endl;
      if (iPivotCol > 0){ // if at least non-zero basis entry found - pivot it
        // Define pivoting multiplier
        double mult = tI[i + 1][iPivotCol];
        //if (mult < -eps){ // make it positive
        //  mult *= -1;
        //}
//        Rcout << "mult=" << mult << std::endl;
        // For the pivot row, divide through the pivot element
        for (int j = 0; j <= tInCol; j++){
          tI[i + 1][j] /= mult;
        }
        // Go through all rows and columns
        for (int j = 0; j <= tInRow; j++){
          // If pivot row - skip, it has been just divided before
          if (j == i + 1){
            continue;
          }
          // Determine row multiplier
          double mult = tI[j][iPivotCol] / tI[i + 1][iPivotCol];
          // Do row operation, i.e. subtract pivot row multiplied with 'mult'
          for (int k = 0; k <= tInCol; k++){
            tI[j][k] -= tI[i + 1][k] * mult;
          }
        }
        // Change the basis record
        tIbasis[i] = iPivotCol;
      }else{
        nZeroBasis++; // one less row will be in the tableau of Phase II
      }
    }
  }
  // 2) Create table structure for Phase II (nrow based on above)
  // Clean the memory if occupied
  if (tII){
    delete[] tII;
    delete[] tIIraw;
  }
//  printTableau(tI, tInRow, tInCol);
//  printBasis(tIbasis, tInRow);
  // Determine size and allocate structures
  tIInRow = tInRow - nZeroBasis;
  tIInCol = AnCol;
  tIIraw = new double[(tIInRow + 1) * (tIInCol + 1)];
  tII = asMatrix(tIIraw, tIInRow + 1, tIInCol + 1); // new
  tIInBasis = tInBasis - nZeroBasis;
  tIIbasis = new int[tIInBasis];
//  Rcout << tIInRow << "x" << tIInCol << "(" << tIInBasis << ")" << std::endl;
  // 3) Copy basis rows (i.e. skip non-basis ones) and basis columns
  //// Copy the goal (zero) row
  //memcpy(tII[0], tI[0], (tIInCol + 1) * sizeof(double));
  // Copy basis rows
  int iRow = 1;
  for (int i = 1; i <= tInRow; i++){ // for each row of the "old" tableau
    if (tIbasis[i - 1] <= AnCol){ // if basis (slack) variable
      // Copy the row and the corresponding basis entry
      memcpy(tII[iRow], tI[i], (tIInCol + 1) * sizeof(double));
      tIIbasis[iRow - 1] = tIbasis[i - 1];
      iRow++;
    }
  }
  // 4) Assemble the goal (zero) row
  for (int i = 0; i <= tIInCol; i++){
//    Rcout << C[i] << " ";
    // Calculate reduced cost of one entry
    double tmp = 0;
    for (int j = 0; j < tIInBasis; j++){
      tmp += C[tIIbasis[j] - 1] * tII[j + 1][i];
    }
    if (i > 0){ // if inner of the tableau
      tII[0][i] = C[i - 1] - tmp;
    }else{
      tII[0][i] = -tmp;
    }
  }
//  Rcout << std::endl;
  return 0;
}

// Read the solution from the class structures
int simplexMethod::readSolution(int *status, double *optVal, double* coefs){
  *status = this->status;
  if (this->status == 0){
    *optVal = this->optVal;
    for (int i = 0; i < AinpnCol; i++){
      coefs[i] = this->coefs[i];
    }
  }
  return 0;
}

// Delete allocated memory
int simplexMethod::freeMemory(){
  if (A){
    delete[] A;
    delete[] Araw;
    delete[] B;
    delete[] C;
    A = 0;
    Araw = 0;
    B = 0;
    C = 0;
  }
  if (Ainp){
    delete[] Ainp;
    delete[] AinpRaw;
    delete[] Binp;
    delete[] Cinp;
    delete[] Dinp;
    Ainp = 0;
    AinpRaw = 0;
    Binp = 0;
    Cinp = 0;
    Dinp = 0;
  }
  if (tI){
    delete[] tI;
    delete[] tIraw;
    delete[] tIbasis;
    tI = 0;
    tIraw = 0;
    tIbasis = 0;
  }
  if (tII){
    delete[] tII;
    delete[] tIIraw;
    delete[] tIIbasis;
    tII = 0;
    tIIraw = 0;
    tIIbasis = 0;
  }
  if (coefs){
    delete[] coefs;
    coefs = 0;
  }
  return 0;
}

int simplexMethod::scale(){
  // Initialize constants
  double abMin = DBL_MAX;
  double abMax = 0;
  // Determine minimum and maximum value
  for (int i = 0; i < AinpnRow; i++){
    // For the constraints matrix
    for (int j = 0; j < AinpnCol; j++){
      if (Ainp[i][j] != 0){ // only non-zero elements participate
        if (fabs(Ainp[i][j]) < abMin){
          abMin = fabs(Ainp[i][j]);
        }
        if (fabs(Ainp[i][j] > abMax)){
          abMax = fabs(Ainp[i][j]);
        }
      }
    }
    // For the right-hand side of the constraints
    if (Binp[i] != 0){
      if (fabs(Binp[i]) < abMin){
        abMin = fabs(Binp[i]);
      }
      if (fabs(Binp[i]) > abMax){
        abMax = fabs(Binp[i]);
      }
    }
  }
  // Determine the scaling constant (and save it)
  abDen = sqrt(abMin * abMax);
//  Rcout << "Scaling constant is " << abDen << std::endl;
  // Scale
  for (int i = 0; i < AinpnRow; i++){
    for (int j = 0; j < AinpnCol; j++){
      Ainp[i][j] /= abDen;
    }
    Binp[i] /= abDen;
  }
  return 0;
}

int simplexMethod::backScale(){
  // Put into the initial scale (working on tII)
  for (int i = 1; i <= tIInRow; i++){
    for (int j = 0; j <= tIInCol; j++){
      tII[i][j] *= abDen;
    }
  }
  return 0;
}
