/*
  File:             radial.cpp
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

void radial(double *yobs, double *xobs, 
            int *m, int *n, int *nobs, double *yref, double *xref, int *nref, 
            int *rts, int *ort, int *ifqh, int *printlevel, double *sol, 
            int *ifsintensities, double *intensities){
//void radial(double *yobs, double *xobs, int *m, int *n, int *nobs, 
//			double *yref, double *xref, int *nref, int *rts, int *ort, int *ifqh, 
//			int *printlevel, double *sol){
 
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
	} else {
		_yref = yref;
		_xref = xref;
		_nref = nref;
	}
	
	/* Create LP-structures */
	int i, j;
//	glp_prob *lp;
//	lp = glp_create_prob();
//	glp_term_out(GLP_OFF);
  simplexMethod sm = simplexMethod();
	/* Add rows */
	int nrow;
	if (*rts == 1 || *rts == 2){
		nrow = *m + *n + 1;
	} else {
		if (*rts == 3){
			nrow = *m + *n;
		}
	}
	int ncol = *_nref + 1;
	int ne = nrow * ncol;
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
	/* Add columns */
//	glp_add_cols(lp, *_nref + 1);
//	glp_set_col_bnds(lp, 0 + 1, GLP_LO, 0., 0.);
//	glp_set_obj_coef(lp, 0 + 1, 1.);
  C[0] = 1;
	for(i = 0; i < *_nref; i++){
//		glp_set_col_bnds(lp, 1 + i + 1, GLP_LO, 0., 0.);
//		glp_set_obj_coef(lp, 1 + i + 1, 0.);
    C[1 + i] = 0;
	}
	/* Set constraint matrix */
	//int *rowIndices = new int[ne];
	//int *colIndices = new int[ne];
	//double *values = new double[ne];
	int *rowIndices = (int *) malloc(ne * sizeof(int));
	int *colIndices = (int *) malloc(ne * sizeof(int));
	double *values = (double *) malloc(ne * sizeof(double));
	for (i = 0; i < nrow; i++){
		for (j = 0; j < ncol; j++){
			rowIndices[i * ncol + j] = i + 1;
			colIndices[i * ncol + j] = j + 1;
			/* Determine the element of the matrix */
			if (i < *m){
				if (j != 0){
//					values[i * ncol + j] = _yref[(j - 1) * *m + i];
          A[i * ncol + j] = _yref[(j - 1) * *m + i];
				}
			}else{
				if (i >= *m && i < *m + *n){
					if (j != 0){
//						values[i * ncol + j] = _xref[(j - 1) * *n + (i - *m)];
            A[i * ncol + j] = _xref[(j - 1) * *n + (i - *m)];
					}
				}else{
					if (i == *m + *n){
						if (j == 0){
//							values[i * ncol + j] = 0.;
              A[i * ncol + j] = 0.;
						}else{
//							values[i * ncol + j] = 1.;
              A[i * ncol + j] = 1.;
						}
					}
				}
			}
		}
	}
	switch (*ort){
	case 1:
		/* Set the direction of optimization */
//		glp_set_obj_dir(lp, GLP_MIN);
		/* Supplement right-hand side */
		for (i = 0; i < *n; i++){
//			glp_set_row_bnds(lp, *m + i + 1, GLP_UP, 0., 0.);
      B[*m + i] = 0;
		  D[*m + i] = 1;
		}
		/* Supplement constraint matrix */
		for (i = 0; i < *m; i++){
//			values[i * ncol] = 0.;
		  A[i * ncol] = 0.;
		}
		break;
	case 2:
		/* Set the direction of optimization */
//		glp_set_obj_dir(lp, GLP_MAX);
    C[0] = -C[0];
		/* Supplement right-hand side */
		for (i = 0; i < *m; i++){
//			glp_set_row_bnds(lp, i + 1, GLP_LO, 0., 0.0);
      B[i] = 0;
		  D[i] = 2;
		}		
		/* Set constraint matrix */
		for (i = *m; i < *m + *n; i++){
//			values[i * ncol] = 0.;
		  A[i * ncol] = 0.;
		}
	}
	/* Run LP-optimization for all observations */
	int index;
	for (index = 0; index < *nobs; index++){
		switch (*ort){
		case 1:
			/* Update right-hand side */
			for (i = 0; i < *m; i++){
//				glp_set_row_bnds(lp, i + 1, GLP_LO, yobs[index * *m + i], 0.);
			  B[i] = yobs[index * *m + i];
			  D[i] = 2;
			}
			/* Update constraint matrix */
			for (i = *m; i < *m + *n; i++){
//				values[i * ncol] = -xobs[index * *n + (i - *m)];
        A[i * ncol] = -xobs[index * *n + (i - *m)];
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
			for (i = 0; i < *m; i++){
//				values[i * ncol] = -yobs[index * *m + i];
        A[i * ncol] = -yobs[index * *m + i];
			}
		}
		/* Load constraint matrix */
//		glp_load_matrix(lp, ne, &rowIndices[-1], &colIndices[-1], &values[-1]);
    sm.inputProblem(A, B, C, D, nrow, ncol);
		// glp_load_matrix(lp, ne, rowIndices, colIndices, values);
		/* Execute linear solver */
//		glp_scale_prob(lp, GLP_SF_GM);
//		glp_simplex(lp, NULL);
		sm.solve();
		int status = -1;
		double optVal = 0;
		double* coefs = new double[ncol];
		sm.readSolution(&status, &optVal, coefs);
		sm.freeMemory();
		if (*ort == 2){
		  optVal = -optVal;
		}
//		if (glp_get_status(lp) == GLP_OPT){
    if (status == 0){
			/* Save solution */
//			sol[index] = glp_get_obj_val(lp);
      sol[index] = optVal;
      if (*ifsintensities){
        // save intensities
        for(i = 0; i < *_nref; i++){
          // begin matter of double precision
          if( fabs ( coefs[i + 1] ) < 1e-8 ){
            intensities[index * *_nref + i] = 0.0;
          } else {
            intensities[index * *_nref + i] = coefs[i + 1];
          }
          // end matter of double precision
        }
      }
		}else{
			/* No solution found */
			sol[index] = 0./0.; // TODO!!!???
		}
		if (*printlevel >= 3){
			Rprintf(" Current observation is %3i/%3i, Farrell measure = %4.4f \n", index + 1, *nobs, sol[index] );
		}
		delete[] coefs;
	}
	/* Release memory */
	// delete[] rowIndices;
	// delete[] colIndices;
	// delete[] values;
	free(rowIndices);
	free(colIndices);
	free(values);	
//	glp_delete_prob(lp);
	if (*ifqh){
		// delete[] _yref;
		// delete[] _xref;
		// delete[] _nref;
		free(_nref);
		free(_yref);
		free(_xref);
	}
	delete[] A;
	delete[] B;
	delete[] C;
	delete[] D;
}

#ifdef __cplusplus
}
#endif
