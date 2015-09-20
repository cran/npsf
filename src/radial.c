#include <R.h>
#include <Rmath.h>
#include "glpk.h"
#include "qhAdapter.h"

void radial(double *yobs, double *xobs, int *m, int *n, int *nobs, 
			double *yref, double *xref, int *nref, int *rts, int *ort, int *ifqh, int *printlevel, 
			double *sol){
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
		convhull(_zref, *nref, *m + *n, vertexIndices);
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
	int i, j;
	glp_prob *lp;
	lp = glp_create_prob();
	glp_term_out(GLP_OFF);
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
	glp_add_rows(lp, nrow);
	if (*rts == 1){
		glp_set_row_bnds(lp, *m + *n + 1, GLP_FX, 1., 1.);
	} else {
		if (*rts == 2){
			glp_set_row_bnds(lp, *m + *n + 1, GLP_UP, 0., 1.);
		}
	}
	/* Add columns */
	glp_add_cols(lp, *_nref + 1);
	glp_set_col_bnds(lp, 0 + 1, GLP_LO, 0., 0.);
	glp_set_obj_coef(lp, 0 + 1, 1.);
	for(i = 0; i < *_nref; i++){
		glp_set_col_bnds(lp, 1 + i + 1, GLP_LO, 0., 0.);
		glp_set_obj_coef(lp, 1 + i + 1, 0.);
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
					values[i * ncol + j] = _yref[(j - 1) * *m + i];
				}
			}else{
				if (i >= *m && i < *m + *n){
					if (j != 0){
						values[i * ncol + j] = _xref[(j - 1) * *n + (i - *m)];
					}
				}else{
					if (i == *m + *n){
						if (j == 0){
							values[i * ncol + j] = 0.;
						}else{
							values[i * ncol + j] = 1.;
						}
					}
				}
			}
		}
	}
	switch (*ort){
	case 1:
		/* Set the direction of optimization */
		glp_set_obj_dir(lp, GLP_MIN);
		/* Supplement right-hand side */
		for (i = 0; i < *n; i++){
			glp_set_row_bnds(lp, *m + i + 1, GLP_UP, 0., 0.);
		}
		/* Supplement constraint matrix */
		for (int i = 0; i < *m; i++){
			values[i * ncol] = 0.;
		}
		break;
	case 2:
		/* Set the direction of optimization */
		glp_set_obj_dir(lp, GLP_MAX);
		/* Supplement right-hand side */
		for (i = 0; i < *m; i++){
			glp_set_row_bnds(lp, i + 1, GLP_LO, 0., 0.0);
		}		
		/* Set constraint matrix */
		for (i = *m; i < *m + *n; i++){
			values[i * ncol] = 0.;
		}
	}
	/* Run LP-optimization for all observations */
	int index;
	for (index = 0; index < *nobs; index++){
		switch (*ort){
		case 1:
			/* Update right-hand side */
			for (i = 0; i < *m; i++){
				glp_set_row_bnds(lp, i + 1, GLP_LO, yobs[index * *m + i], 0.);
			}
			/* Update constraint matrix */
			for (i = *m; i < *m + *n; i++){
				values[i * ncol] = -xobs[index * *n + (i - *m)];
			}
			break;
		case 2:
			/* Update right-hand side */
			for (i = 0; i < *n; i++){
				glp_set_row_bnds(lp, *m + i + 1, GLP_UP, 0., xobs[index * *n + i]);
			}
			/* Update constraint matrix */
			for (i = 0; i < *m; i++){
				values[i * ncol] = -yobs[index * *m + i];
			}
		}
		/* Load constraint matrix */
		glp_load_matrix(lp, ne, &rowIndices[-1], &colIndices[-1], &values[-1]);
		// glp_load_matrix(lp, ne, rowIndices, colIndices, values);
		/* Execute linear solver */
		glp_scale_prob(lp, GLP_SF_GM);
		glp_simplex(lp, NULL);
		if (glp_get_status(lp) == GLP_OPT){
			/* Save solution */
			sol[index] = glp_get_obj_val(lp);
		}else{
			/* No solution found */
			sol[index] = 0./0.;
		}
		if (*printlevel >= 3){
			Rprintf(" Current observation is %3i/%3i, Farrell measure = %4.4f \n", index + 1, *nobs, sol[index] );
		}
	}
	/* Release memory */
	// delete[] rowIndices;
	// delete[] colIndices;
	// delete[] values;
	free(rowIndices);
	free(colIndices);
	free(values);	
	glp_delete_prob(lp);
	if (*ifqh){
		// delete[] _yref;
		// delete[] _xref;
		// delete[] _nref;
		free(_nref);
		free(_yref);
		free(_xref);
	}
}
