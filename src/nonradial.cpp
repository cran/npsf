/*
 File:             nonradial.cpp
 Created by:       Pavlo Mozharovskyi, Oleg Badunenko
 First published:  2015-09-20
 Last revised:     2020-06-24
 
 Solves LPs
 */

#include <R.h>
#include <Rmath.h>
//#include "glpk.h"
//#include "qhAdapter.h"
#include "simplexMethod.h"

#define EOF (-1)

#ifdef __cplusplus
extern "C" {
#endif
 
void nonradial(double *yobs, double *xobs, int *m, int *n, int *nobs, 
                double *yref, double *xref, int *nref,  
                int *rts, int *ort, int *ifqh, int *printlevel, 
                double *sol, double *solfull, double *intensities, 
                int *ifsolfull, int *lmdConstr ){
  //void nonradial(double *yobs, double *xobs, int *m, int *n, int *nobs,
  //			double *yref, double *xref, int *nref, int *rts, int *ort, int *ifqh, 
  //      int *printlevel, double *sol){
  
  /*
   // Begin Scaling
   
   // Total columns to scale
   int numColsToScale = 2 * *m + 2 * *n;
   // Columns' pointers
   double **ptrsCol = (double **) malloc(numColsToScale * sizeof(double *));
   // Columns' lengths
   int *numsRow = (int *) malloc(numColsToScale * sizeof(int));
   // Columns' numbers
   int *numsCol = (int *) malloc(numColsToScale * sizeof(int));
   // Create objects with addresses of rows and columns of XOBS etc.
   // Outputs
   for (int i = 0; i < *m; i++){
   // i-th column of yobs
   ptrsCol[i * 2] = yobs + i;
   numsRow[i * 2] = *nobs;
   numsCol[i * 2] = *m;
   // i-th column of yref
   ptrsCol[i * 2 + 1] = yref + i;
   numsRow[i * 2 + 1] = *nref;
   numsCol[i * 2 + 1] = *m;
   }
   // Inputs
   for (int i = 0; i < *n; i++){
   // i-th column of xobs
   ptrsCol[2 * *m + i * 2] = xobs + i;
   numsRow[2 * *m + i * 2] = *nobs;
   numsCol[2 * *m + i * 2] = *n;
   // i-th column of xref
   ptrsCol[2 * *m + i * 2 + 1] = xref + i;
   numsRow[2 * *m + i * 2 + 1] = *nref;
   numsCol[2 * *m + i * 2 + 1] = *n;
   }
   // Namely scaling
   for (int i = 0; i < numColsToScale; i++){
   // Calculate mean
   double tmpMean = 0;
   for (int j = 0; j < numsRow[i]; j++){
   tmpMean += *(ptrsCol[i] + j * numsCol[i]);
   }
   tmpMean /= numsRow[i];
   // Scale if mean significanly larger zero
   if (tmpMean > 1e-8){
   for (int j = 0; j < numsRow[i]; j++){
   *(ptrsCol[i] + j * numsCol[i]) /= tmpMean;
   }
   }
   }
   // Clean the memory
   free(ptrsCol);
   free(numsRow);
   free(numsCol);
   
   // End Scaling
   */
  
  /* Temporary structures for reference data */
  double *_yref;
  double *_xref;
  int *_nref;
  if (*ifqh){
   /* Filter reference data if necessary s.t. only convex hull leaves */
   //double *_zref = new double[(*m + *n) * *nref];
   double *_zref = (double *) malloc(((*m + *n) * *nref) * sizeof(double));
   for (int i = 0; i < *nref; i++){
    for (int j = 0; j < (*m + *n); j++){
     if (j < *m){
      _zref[i * (*m + *n) + j] = yref[i * *m + j];
     }else{
      _zref[i * (*m + *n) + j] = xref[i * *n + (j - *m)];
     }
    }
   }
   //int *vertexIndices = new int[*nref + 1];
   int *vertexIndices = (int *) malloc((*nref + 1) * sizeof(int));
   //convhull(_zref, *nref, *m + *n, vertexIndices);
   //CH_qhull(_zref, *nref, *m + *n, vertexIndices);
   _nref = vertexIndices;
   free(_zref);
   //_yref = new double[*m * *_nref];
   //_xref = new double[*n * *_nref];
   _yref = (double *) malloc((*m * *_nref) * sizeof(double));
   _xref = (double *) malloc((*n * *_nref) * sizeof(double));
   for (int i = 0; i < vertexIndices[0]; i++){
    for (int j = 0; j < *m; j++){
     _yref[i * *m + j] = yref[vertexIndices[i + 1] * *m + j];
    }
    for (int j = 0; j < *n; j++){
     _xref[i * *n + j] = xref[vertexIndices[i + 1] * *n + j];
    }
   }
  }else{
   _yref = yref;
   _xref = xref;
   _nref = nref;
  }
  /* Create LP-structures */
  // Rprintf("this is 0\n");
  int i, j;
  //	glp_prob *lp;
  //	lp = glp_create_prob();
  //	glp_term_out(GLP_OFF);
  simplexMethod sm = simplexMethod();
  /* Add rows */
  int nrow, ncol1, yxRTSlength;
  if (*rts == 1 || *rts == 2){
   nrow = *m + *n + 1;
  }else{
   if (*rts == 3){
    nrow = *m + *n;
   }
  }
  yxRTSlength = nrow;
  // Rprintf("this is 1\n");
  // int ncol = *_nref + 1;
  switch (*ort){
  case 1:
   ncol1 = *n;
   break;
  case 2:
   ncol1 = *m;
  }
  
  // begin adding rows if lmdConstr
  if(*lmdConstr){
    nrow += ncol1;
  }
  // now {nrow = *m + *n + ncol1}
  // end adding rows if lmdConstr
  
  // double ncol2 = ncol1;
  // Rprintf("int = %2i, double = %4.2f\n", ncol1, ncol2);
  // Rprintf("nrow = %2i, ncol = %2i\n", nrow, ncol1);
  int ncol = ncol1 + *_nref;
  int ne = nrow * ncol;
  // Rprintf("nrow = %2i, ncol = %2i, ncol1 = %2i, *_nref = %2i, ne = %2i\n", 
          // nrow, ncol, ncol1, *_nref, ne);
  // Rprintf("yxRTSlength = %2i, *lmdConstr = %2i, *ifqh = %2i, *ifsolfull = %2i, *printlevel = %2i\n", 
          // yxRTSlength, *lmdConstr, *ifqh, *ifsolfull, *printlevel);
  //	glp_add_rows(lp, nrow);
  double* A = new double[ne]; // left-hand-side matrix of the constraints
  double* B = new double[nrow]; // right-hand side of the constraints
  double* C = new double[ncol]; // the goal function
  int* D = new int[nrow]; // (in)equality directions
  if (*rts == 1){
   //		glp_set_row_bnds(lp, *m + *n + 1, GLP_FX, 1., 1.);
   B[*m + *n] = 1;
   D[*m + *n] = 0;
  } else {
   if (*rts == 2){
    //			glp_set_row_bnds(lp, *m + *n + 1, GLP_UP, 0., 1.);
    B[*m + *n] = 1;
    D[*m + *n] = 1;
   }
  }
  // Rprintf("this is 2\n");
  /* Add columns */
  //	glp_add_cols(lp, ncol);
  // second part of the objective f-n
  for(i = 0; i < *_nref; i++){
   //		glp_set_col_bnds(lp, ncol1 + i + 1, GLP_LO, 0., 0.);
   //		glp_set_obj_coef(lp, ncol1 + i + 1, 0.);
   C[ncol1 + i] = 0;
  }
  // Rprintf("this is 3\n");
  /* Set constraint matrix */
  //int *rowIndices = new int[ne];
  //int *colIndices = new int[ne];
  //double *values = new double[ne];
  int *rowIndices = (int *) malloc(ne * sizeof(int));
  int *colIndices = (int *) malloc(ne * sizeof(int));
  double *values = (double *) malloc(ne * sizeof(double));
  // Rprintf("this is 4\n");
  // make them all zeros first and fill in indices
  for (i = 0; i < nrow; i++){
   for (j = 0; j < ncol; j++){
    rowIndices[i * ncol + j] = i + 1;
    colIndices[i * ncol + j] = j + 1;
    //			values[i * ncol + j] = 0.0;
    A[i * ncol + j] = 0;
   }
  }
  // fill in the second/right part of the constraint matrix
  // @@ instead of {i < nrow}, we now use {i < yxRTSlength}
  // @@ because we do not touch lower panel
  for (i = 0; i < yxRTSlength; i++){
   for (j = ncol1 - 1; j < ncol; j++){
    // rowIndices[i * ncol + j] = i + 1;
    // colIndices[i * ncol + j] = j + 1;
    /* Determine the element of the matrix */
    if (i < *m){
     if (j > ncol1-1){
      //					values[i * ncol + j] = _yref[(j - ncol1) * *m + i];
      A[i * ncol + j] = _yref[(j - ncol1) * *m + i];
     }
    }else{
     if (i >= *m && i < *m + *n){
      if (j > ncol1-1){
       //						values[i * ncol + j] = _xref[(j - ncol1) * *n + (i - *m)];
       A[i * ncol + j] = _xref[(j - ncol1) * *n + (i - *m)];
      }
     }else{
      if (i == *m + *n){
       if (j <= ncol1-1){
        //							values[i * ncol + j] = 0.;
        A[i * ncol + j] = 0;
       }else{
        //							values[i * ncol + j] = 1.;
        A[i * ncol + j] = 1.;
       }
      }
     }
    }
   }
  }
  // Rprintf("this is 4\n");
  // @@ this is where we touch lower panel
  // begin adding 1 or -1 on diagonal of size ncol1 if lmdConstr; B; D
  if(*lmdConstr){
    
    // begin A
    // in input base A[.] <= 1
    // in output base -A[.] <= -1
    // switch (*ort){
    // case 1:
    //   for (i = 0; i < ncol1; i++){
    //     A[ (i + yxRTSlength) * ncol + i] = 1.0;
    //     // Rprintf("i = %2i,(i+yxRTSlength) * ncol + i = %2i\n", i, (i+yxRTSlength) * ncol + i);
    //   }
    //   break;
    // case 2:
    //   for (i = 0; i < ncol1; i++){
    //     A[ (i + yxRTSlength) * ncol + i] = 1.0;
    //     // Rprintf("i = %2i, (i+yxRTSlength) * ncol + i = %2i\n", i, (i+yxRTSlength) * ncol + i);
    //   }
    // }
    for (i = 0; i < ncol1; i++){
      A[ (i + yxRTSlength) * ncol + i] = 1.0;
      // Rprintf("i = %2i,(i+yxRTSlength) * ncol + i = %2i\n", i, (i+yxRTSlength) * ncol + i);
    }
    // end A
    
    // either greater than or less than or equal 1 for all lambdas
    // begin B
    // in input base A[.] <= 1
    // in output base -A[.] <= -1
    // switch (*ort){
    // case 1:
    //   for (i = yxRTSlength; i < nrow; i++){
    //     B[i] = 1.0;
    //   }
    //   break;
    // case 2:
    //   for (i = yxRTSlength; i < nrow; i++){
    //     B[i] = 1.0;
    //   }
    // }
    
    for (i = yxRTSlength; i < nrow; i++){
      B[i] = 1.0;
    }
    
    // end B
    // in input base A[.] <= 1
    // in output base -A[.] <= -1
    // so always less or equal
    
    // begin D
    // for (i = yxRTSlength; i < nrow; i++){
    //   D[i] = 1;
    // }
    switch (*ort){
    case 1:
      for (i = yxRTSlength; i < nrow; i++){
        D[i] = 1;
      }
      break;
    case 2:
      for (i = yxRTSlength; i < nrow; i++){
        D[i] = 2;
      }
    }
    // end D
  }
  // end adding 1 or -1 on diagonal of size ncol1 if lmdConstr; B; D
  
  // Rprintf("this is 5\n");
  switch (*ort){
  case 1:
   /* Set the direction of optimization */
   //		glp_set_obj_dir(lp, GLP_MIN);
   /* Supplement right-hand side */
   for (i = 0; i < *n; i++){
    //			glp_set_row_bnds(lp, *m + i + 1, GLP_UP, 0.0, 0.0);
    B[*m + i] = 0;
    D[*m + i] = 1;
   }
   /* Supplement constraint matrix */
   // for (int i = 0; i < *m; i++){
   // 	values[i * ncol] = 0.;
   // }
   break;
  case 2:
   /* Set the direction of optimization */
   //		glp_set_obj_dir(lp, GLP_MAX);
   /* Supplement right-hand side */
   for (i = 0; i < *m; i++){
    //			glp_set_row_bnds(lp, i + 1, GLP_LO, 0.0, 0.0);
    B[i] = 0;
    D[i] = 2;
   }
   /* Set constraint matrix */
   // for (int i = *m; i < *m + *n; i++){
   // 	values[i * ncol] = 0.;
   // }
  }
  // Rprintf("this is 6\n");
  /* Run LP-optimization for all observations */
  int index;
  double zeros;
  int reversed = 0;
  for (index = 0; index < *nobs; index++){
   zeros = 0.0;
   // if reversed for previous i, get the signs right (begin)
   // no need to check lmdContr as 'reversed = 1 can only be if lmdContr
   if(reversed == 1){
     switch (*ort){
     case 1:
       for (i = yxRTSlength; i < nrow; i++){
         D[i] = 1;
       }
       break;
     case 2:
       for (i = yxRTSlength; i < nrow; i++){
         D[i] = 2;
       }
     }
   }
   // if reversed for previous i, get the signs right (end)
   switch (*ort){
   case 1:
    /* Update right-hand side */
    for (i = 0; i < *m; i++){
     //				glp_set_row_bnds(lp, i + 1, GLP_LO, yobs[index * *m + i], 0.);
     B[i] = yobs[index * *m + i];
     D[i] = 2;
    }
    /* Update constraint matrix */
    // zero in objective f-n if x == 0
    // count zeros
    for(i = 0; i < ncol1; i++){
     //				glp_set_col_bnds(lp, i + 1, GLP_LO, 0., 0.);
     if(xobs[index * *n + (i - 0)] == 0){
      zeros += 1.0;
      //				glp_set_obj_coef(lp, i + 1, 0.);
      C[i] = 0;
     } else {
      //					glp_set_obj_coef(lp, i + 1, 1.);
      C[i] = 1.;
     }
    }
    // Rprintf("index=%2i, nu of zeros=%2i \n", index, zeros);
    // diagonal elements of (*n)x(*n) matrix
    for (i = *m; i < *m + *n; i++){
     //				values[i * ncol + (i - *m)] = -xobs[index * *n + (i - *m)];
     A[i * ncol + (i - *m)] = -xobs[index * *n + (i - *m)];
      // if (index == 0){
      //  Rprintf("i = %2i,i * ncol + (i - *m) = %2i\n", i, i * ncol + (i - *m));
      // }
    }
    break;
   case 2:
    /* Update right-hand side */
    for (i = 0; i < *n; i++){
     //				glp_set_row_bnds(lp, *m + i + 1, GLP_UP, 0., xobs[index * *n + i]);
     B[*m + i] = xobs[index * *n + i];
     D[*m + i] = 1;
    }
    /* Update constraint matrix */
    // zero in objective f-n if x == 0
    // count zeros
    for(i = 0; i < *m; i++){
     //				glp_set_col_bnds(lp, i + 1, GLP_LO, 0., 0.);
     if(yobs[index * *m + i] == 0){
      zeros += 1.0;
      //					glp_set_obj_coef(lp, i + 1, 0.);
      C[i] = 0;
     } else {
      //					glp_set_obj_coef(lp, i + 1, 1.);
      C[i] = -1.;
     }
     // diagonal elements of (*n)x(*n) matrix
     //				values[i * ncol + i] = -yobs[index * *m + i];
     A[i * ncol + i] = -yobs[index * *m + i];
     //if (index == 0){
       // Rprintf("i = %2i,i * ncol + i = %2i\n", i, i * ncol + i);
     //}
    }
   }
   /* Load constraint matrix */
   // glp_load_matrix(lp, ne, rowIndices, colIndices, values);
   //		glp_load_matrix(lp, ne, &rowIndices[-1], &colIndices[-1], &values[-1]);
   // if(index == 0){
   //   for (i = 0; i < ncol; i++){
   //     Rprintf("i = %2i,C[i] = %4.2f\n", i, C[i]);
   //   }
   // }
   sm.inputProblem(A, B, C, D, nrow, ncol);
   //* Execute linear solver */
   //		glp_scale_prob(lp, GLP_SF_GM);
   //		glp_simplex(lp, NULL);
   // sm.solve();
   
   
   // prepare for status and solution reading
   int status = -1;
   double optVal = 0;
   double* coefs = new double[ncol];
   
   /* In the unconstrained case, perfom in one step */
   if (*lmdConstr == 0){
     //* Execute linear solver */
     sm.solve();
     sm.readSolution(&status, &optVal, coefs);
   } else {
     /* In the constrained case proceed in two steps */
     /* First, try to solve pretending that the point is not outside the frontier */

     // solve
     sm.solve();
     
     // read solution and status
     sm.readSolution(&status, &optVal, coefs);
     // Rprintf("Status for i = %2i, is %2i\n", index+1, status);
     
     if ( status != 0 ){
       reversed = 1;
       /* Second, reverse the constraints around one */
       // reverse D (begin)
       // Rcout << "Had to reverse for" << index << "";
       // Rprintf("Had to reverse for i = %2i\n", index+1);
       // Rcout << Dinp[i] << " ";
       // Rcout << std::endl;
       switch (*ort){
       case 1:
         for (i = yxRTSlength; i < nrow; i++){
           D[i] = 2;
         }
         break;
       case 2:
         for (i = yxRTSlength; i < nrow; i++){
           D[i] = 1;
         }
       }
       // reverse D (end)
       sm.inputProblem(A, B, C, D, nrow, ncol);
       sm.solve();
       sm.readSolution(&status, &optVal, coefs);
       // Rprintf("Status for i = %2i, is %2i\n", index+1, status);
     }
   }
   
   
  
   // int status = -1;
   // double optVal = 0;
   // double* coefs = new double[ncol];
   // sm.readSolution(&status, &optVal, coefs);
   //sm.freeMemory();
   if (*ort == 2){
    optVal = -optVal;
   }
   // Rprintf("i = %2i, status = %2i, optVal = %4.4f\n", index+1, status, optVal/ncol1);
   //		if (glp_get_status(lp) == GLP_OPT){
   if (status == 0){
    /* Save solution */
    //			sol[index] = (glp_get_obj_val(lp) + zeros) / (double) ncol1;
    sol[index] = (optVal + zeros) / (double) ncol1;
     if (*ifsolfull){
     /* Save detailed solution */
     switch (*ort){
     case 1:
      /* Consider zeros */
      for(i = 0; i < *n; i++){
       if(xobs[index * *n + (i - 0)] == 0){
        solfull[index * *n + (i - 0)] = 0.0;
       } else {
        solfull[index * *n + (i - 0)] = coefs[i];
        // begin matter of double precision
        if( fabs (1 - coefs[i] ) < 1e-8 ){
         solfull[index * *n + (i - 0)] = 1.0;
        }
        // end matter of double precision
       }
      }
      break;
     case 2:
      /* Consider zeros */
      for(i = 0; i < *m; i++){
       if(yobs[index * *m + i] == 0){
        solfull[index * *m + i] = 0.0;
       } else {
         solfull[index * *m + i] = coefs[i];
         // begin matter of double precision
         if( fabs (1 - coefs[i] ) < 1e-8 ){
           solfull[index * *m + i] = 1.0;
         }
         // end matter of double precision
       }
      }
     }
    }
     // save intensities
     for(i = 0; i < *_nref; i++){
       // begin matter of double precision
       if( fabs ( coefs[i + ncol1] ) < 1e-8 ){
         intensities[index * *_nref + i] = 0.0;
       } else {
         intensities[index * *_nref + i] = coefs[i + ncol1];
       }
       // end matter of double precision
     }
   }else{
    /* No solution found */
    sol[index] = -999; // TODO!!!???
   }
   if (*printlevel >= 3){
    Rprintf(" Current observation is %3i/%3i, Russell measure = %4.4f \n", index + 1, *nobs, sol[index] );
   }
   delete[] coefs;
  }
  /* Release memory */
  //delete[] rowIndices;
  //delete[] colIndices;
  //delete[] values;
  free(rowIndices);
  free(colIndices);
  free(values);
  //	glp_delete_prob(lp);
  if (*ifqh){
   free(_yref);
   free(_xref);
   free(_nref);
  }
  delete[] A;
  delete[] B;
  delete[] C;
  delete[] D;
 }
 
#ifdef __cplusplus
}
#endif
