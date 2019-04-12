/*
 File:             fourcomp.c
 Created by:       Oleg Badunenko
 First published:  2019-04-12
 Last revised:     2019-04-12
 
 cacluates ll, grad, and hessian of the homoskedastic 4 components model
 */

#include <R.h>
#include <Rmath.h>

// fun0:        Product by ids: NT in; N out
// NT:          zit
// N:           ids (identifiers of IDs)
// 1:           NT, N
// zi:      N

void it_prods(double *zit,
double *idvar, int *NT,
double *ids, int *N,
double *zi){
	for(int i = 0; i < *N; i++) // over panel id, N
	{
		// Need within id; for this we need
		// one more cicle over NT, inside of circle over NT
		zi[i] = 1;
		for(int j = 0; j < *NT; j++)
		{
			if(idvar[j] == ids[i])
			{
				zi[i] *= zit[j];
			}
		}
	}
}

void it_sums(double *zit,
double *idvar, int *NT,
double *ids, int *N,
double *zi){
	for(int i = 0; i < *N; i++) // over panel id, N
	{
		// Need within id; for this we need
		// one more cicle over NT, inside of circle over NT
		zi[i] = 0;
		for(int j = 0; j < *NT; j++)
		{
			if(idvar[j] == ids[i])
			{
				zi[i] += zit[j];
			}
		}
	}
}


// in:
// NxR:	W
// NTxK:	Z
// NT:	C, idvar
// N:		ids,
// K4:	theta = c(gam, lmd, sig, sigw, sigh)
// 1:		NT, N, R, idmax, 

// out
// 1:		lnls
// K4:	grad
// (( k4^2 + k4 ) / 2):	hess (lower triangle)

// prod == 1  => production function
// prod == -1 => cost function

void gtre(int *prod, double *W, double *H, int *N, int *R, double *Z, int *NT, double *C, 
double *ids, double *idvar, int *idlenmax, double *theta, int *K4, int *BHHH,
double *lnls, double *grad, double *hesstri){
	int k4H = ( pow(K4[0],2) + K4[0])/2;
	// calculate CZg = C - Z*g
	// Allocate memory
	double * CZg   = (double *) malloc (*NT * sizeof(double));
	// Multiplication of Z and theta
	for(int q = 0; q < *NT; q++){
		for(int w = 0; w < *K4-4; w++){
			if(w == 0){
				CZg[q] = C[q];
			}
			CZg[q] -= Z[ w*NT[0] + q ] * theta[w];
		}
		// Rprintf("q=%2i, CZg=%4.3E \n", q, CZg[q]);
	}
	// Allocate memory
	double * Pir    = (double *) malloc (*R  * sizeof(double));
	double * Pi     = (double *) malloc (*N  * sizeof(double));
	double * gitr   = (double *) malloc (*K4 * sizeof(double));
	double * gir    = (double *) malloc (R[0]*K4[0] * sizeof(double));
	double * gibar  = (double *) malloc (N[0]*K4[0] * sizeof(double));
	double * Hitr   = (double *) malloc (k4H * sizeof(double));
	double * Hitr2  = (double *) malloc (k4H * sizeof(double));
	double * Hitr3  = (double *) malloc (k4H * sizeof(double));
	double * Hir    = (double *) malloc (R[0]*k4H * sizeof(double));
	double * Hir2   = (double *) malloc (R[0]*k4H * sizeof(double));
	// double * Hir3   = (double *) malloc (R[0]*k3H * sizeof(double));
	double * Hibar  = (double *) malloc (N[0]*k4H * sizeof(double));
	double * Hibar2 = (double *) malloc (N[0]*k4H * sizeof(double));
	// double * Hibar3 = (double *) malloc (N[0]*k3H * sizeof(double));
	double * HibarBHHH = (double *) malloc (k4H * sizeof(double));
		
	
	
	int k, r, ijmiss, i, j, ii, jj;
	// get the sums of Pir first to get Qir later
	for(i = 0; i < *N; i++){
		// now within r
		// initialize for likelihood
		Pi[i] = 0.0;
		// go on
		for(r = 0; r < *R; r++){
			// now within i
			// int Ti = 0;
			// initialize for likelihood
			Pir[r] = 1.0;
			// go on
			for(j = 0; j < *NT; j++)
			{
				if(idvar[j] == ids[i])
				{
					// here is the main action
					double A = (CZg[j] - theta[*K4-2] * W[r*N[0] + i ] + theta[*K4-1] * H[r*N[0] + i ]*prod[0] ) / theta[*K4-3];
					double B = A * theta[*K4-4];
					// begin calculate "itr" element
					double Pitr = 2.0 * pow(theta[*K4-3],-1.0) * dnorm(A, 0, 1, 0) * pnorm(-B*prod[0], 0, 1, 1, 0);
					// end calculate "itr" element
					// product over "t"
					Pir[r] *= Pitr;
				} // end condition "idvar[j] == ids[i]"
			} // end checking if in this panel
			// sum of Pir
			Pi[i] += Pir[r];
		} // end of loop within "r"
		// Rprintf("i=%2i, Pi=%4.3E \n", i, Pi[i]);
	} // end of loop within "i"
	
	// Now the code for calculation of lnls, gr and hess
	
	// initialise those that come from R
	*lnls = 0.0;
	
	for(k = 0; k < *K4; k++){
		grad[k] = 0.0;
	}
	for(k = 0; k < k4H; k++){
		hesstri[k] = 0.0;
	}
	for(i = 0; i < *N; i++){
		// now within r
		// initialize for gradient
		for(k = 0; k < *K4; k++){
			gibar[k*N[0] + i] = 0.0;
		}
		if(*BHHH == 0){
			for(k = 0; k < k4H; k++){
				Hibar[k*N[0] + i] = 0.0;
				Hibar2[k*N[0] + i] = 0.0;
				// Hibar3[k*N[0] + i] = 0.0;
			}
		} // end if_bhhh
		// go on
		for(r = 0; r < *R; r++){
			// now within i
			// int Ti = 0;
			// initialize for likelihood
			Pir[r] = 1.0;
			// initialize for gradient
			for(k = 0; k < *K4; k++){
				gir[k*R[0] + r] = 0.0;	
			}
			if(*BHHH == 0){
				// initialize for hessian
				for(k = 0; k < k4H; k++){
					Hir[k*R[0] + r] = 0.0;
					Hir2[k*R[0] + r] = 0.0;
					// Hir3[k*R[0] + r] = 0.0;
				}
			} // end if_bhhh
			// go on
			for(j = 0; j < *NT; j++)
			{
				if(idvar[j] == ids[i])
				{
					// here is the main action
					double A = (CZg[j] - theta[*K4-2] * W[r*N[0] + i ] + theta[*K4-1] * H[r*N[0] + i ]*prod[0] ) / theta[*K4-3];
					double B = A * theta[*K4-4];
					double C = exp( dnorm(-B*prod[0], 0, 1, 1) - pnorm(-B*prod[0], 0, 1, 1, 1));
					double D = B - C*prod[0];
					double E = pow(theta[*K4-3],-1.0) * (A + theta[*K4-4]*C*prod[0]);
					double F = -1.0 * pow(theta[*K4-3],-2.0) * (1 - pow(theta[*K4-4], 2.0) * C *prod[0] * D);
					double G = C *prod[0] * (1.0 - B * D);
					// @@@ likelihood function
					// begin calculate "itr" element
					double Pitr = 2.0 * pow(theta[*K4-3],-1.0) * dnorm(A, 0, 1, 0) * pnorm(-B*prod[0], 0, 1, 1, 0);
					// end calculate "itr" element
					// Ti++;
					// if(i == 1){
					// 	Rprintf("r = %2i, j = %2i, i = %2i, Pitr = %4.8f \n", r, j, i, Pi);
					// }
					// product over "t"
					Pir[r] *= Pitr;
					// @@@ gradient
					// begin calculate "itr" element
					// w.r.t. gammas
					for(k = 0; k < *K4-4; k++){
						gitr[k] = E * Z[ k*NT[0] + j ];
					}
					// w.r.t. lambda
					gitr[*K4-4] = -C *prod[0] * A;
					// w.r.t. sigma
					gitr[*K4-3] = -1.0 * pow(theta[*K4-3],-1.0) + E * A;
					// w.r.t. sigma_w
					gitr[*K4-2] = E * W[r*N[0] + i ];
					// w.r.t. sigma_h
					gitr[*K4-1] = E * -H[r*N[0] + i ]*prod[0];
					// end calculate "itr" element
					// sum over "t"
					for(k = 0; k < *K4; k++){
						gir[k*R[0] + r] += gitr[k];	
					}
					// @@@ hessian
					if(*BHHH == 0){
						// begin first term of the hessian
						// begin calculate "itr" element
						ijmiss = 0.0;
						for(jj = 0; jj < *K4-4; jj++){
							for(ii = 0; ii < *K4-4; ii++){
								if(ii >= jj){
									// begin fill in the ZZ'
									Hitr[ jj * K4[0] + ii - ijmiss] = F * Z[ ii*NT[0] + j ] * Z[ jj*NT[0] + j ];
									// end fill in the ZZ'
									// w.r.t gamma and lambda: sig1^-1 * G1[i,r] * zit[i,]
									Hitr[ jj * K4[0] + ii + 1 - ijmiss] = pow(theta[*K4-3],-1.0) * G * Z[ jj*NT[0] + j];
									// w.r.t gamma and sigma: -1 * sig1^-2 * (2*A1[i,r] + lmd1 * G1[i,r]) * zit[i,]
									Hitr[ jj * K4[0] + ii + 2 - ijmiss] = -1.0 * pow(theta[*K4-3],-2.0) * (2.0*A + theta[*K4-4]*G) * Z[jj*NT[0]+j];
									// w.r.t gamma and sigmaw: F1[i,r] * Wexpanded[i,r] * zit[i,]
									Hitr[ jj * K4[0] + ii + 3 - ijmiss] = F * W[r*N[0] + i ] * Z[jj*NT[0]+j];
									// w.r.t gamma and sigmaw: F1[i,r] * -Hexpanded[i,r] * zit[i,]
									Hitr[ jj * K4[0] + ii + 4 - ijmiss] = F * -H[r*N[0] + i ]*prod[0] * Z[jj*NT[0]+j];
								} else {
									ijmiss++;
								}
							} //end "ii"
						}//end "jj"
						// w.r.t. lambda and lambda: A1[i,r]^2 * C1[i,r] * D1[i,r]
						Hitr[ k4H - 10] = pow(A, 2.0) * C *prod[0] * D;
						// w.r.t. lambda and sigma: A1[i,r] * sig1^-1 * G1[i,r]
						Hitr[ k4H - 9] = A * pow(theta[*K4-3],-1.0) * G;
						// w.r.t lambda and sigmaw: sig1^-1 * Wexpanded[i,r] * G1[i,r]
						Hitr[ k4H - 8] = pow(theta[*K4-3],-1.0) * W[r*N[0] + i ] * G;
						// w.r.t lambda and sigmah: sig1^-1 * -Hexpanded[i,r] * G1[i,r]
						Hitr[ k4H - 7] = pow(theta[*K4-3],-1.0) * -H[r*N[0] + i ]*prod[0] * G;
						// w.r.t. sigma and sigma: sig1^-2 * (1 - 3*A1[i,r]^2 - B1[i,r]*C1[i,r] - B1[i,r]*G1[i,r])
						Hitr[ k4H - 6] = pow(theta[*K4-3],-2.0) * (1 - 3 * pow(A, 2.0) - B*C*prod[0] - B*G);
						// w.r.t gamma and sigmaw: -1 * sig1^-2 * (2*A1[i,r] + lmd1 * G1[i,r]) * Wexpanded[i,r]
						Hitr[ k4H - 5] = -1.0 * pow(theta[*K4-3],-2.0) * (2.0*A + theta[*K4-4]*G) * W[r*N[0] + i ];
						// w.r.t gamma and sigmah: -1 * sig1^-2 * (2*A1[i,r] + lmd1 * G1[i,r]) * -Hexpanded[i,r]
						Hitr[ k4H - 4] = -1.0 * pow(theta[*K4-3],-2.0) * (2.0*A + theta[*K4-4]*G) * -H[r*N[0] + i ]*prod[0];
						// w.r.t sigmaw and sigmaw: F1[i,r] * Wexpanded[i,r]^2
						Hitr[ k4H - 3] = F * pow(W[r*N[0] + i ], 2.0);
						// w.r.t sigmaw and sigmah: -F1[i,r] * Wexpanded[i,r] * Hexpanded[i,r]
						Hitr[ k4H - 2] = -F * W[r*N[0] + i ] * H[r*N[0] + i ]*prod[0];
						// w.r.t sigmah and sigmah: F1[i,r] * Hexpanded[i,r]^2
						Hitr[ k4H - 1] = F * pow(H[r*N[0] + i ], 2.0);
						// if(i <= 1){
						// 	Rprintf("i = %2i, H_sig_u0 = %4.8e \n", i, Hitr[ k4H - 1]);
						// }
						// end calculate "itr" element
						// sum over "t"
						for(k = 0; k < k4H; k++){
							Hir[k*R[0] + r] += Hitr[k];
						}
						// end first term of the hessian
					} // end if_bhhh
				} // end condition "idvar[j] == ids[i]"
			} // end checking if in this panel
			// weighted sums for gradient, sums of columns of (supposedly matrix) gir
			for(k = 0; k < *K4; k++){
				gibar[k*N[0] + i] += gir[k*R[0] + r] * Pir[r] / Pi[i];
			}
			if(*BHHH == 0){
				// weighted sums for hessian, sums of columns of (supposedly matrix) Hir
				for(k = 0; k < k4H; k++){
					Hibar[k*N[0] + i] += Hir[k*R[0] + r] * Pir[r] / Pi[i];
					// Hibar2[k*N[0] + i] += Hir2[k*R[0] + r];
					// Hibar3[k*N[0] + i] += Hir2[k*R[0] + r];
				}
				// begin second and third terms of the hessian
				// begin calculate "ir" element
				ijmiss = 0.0;
				for(jj = 0; jj < *K4; jj++){
					for(ii = 0; ii < *K4; ii++){
						if(ii >= jj){
							Hitr2[ jj * K4[0] + ii - ijmiss] = gir[ii*R[0] + r] * gir[jj*R[0] + r] * Pir[r] / Pi[i];
							// Hitr3[ jj * K3[0] + ii - ijmiss] = gir[ii*R[0] + r] * gir[jj*R[0] + r] * Pir[r] / Pi[i] * Pir[r] / Pi[i];
							// if(i <= 1){
							// 	Rprintf("i = %2i, Hibar2 = %4.8e, Hibar3 = %4.8e \n", i, Hitr2[ jj * K3[0] + ii - ijmiss], Hitr3[ jj * K3[0] + ii - ijmiss]);
							// }
						} else {
							ijmiss++;
						}
					} //end "ii"
				}//end "jj"
				// end calculate "ir" element
			} // end if_bhhh
			// sum over "r"
			if(*BHHH == 0){
				for(k = 0; k < k4H; k++){
					Hibar2[k*N[0] + i] += Hitr2[k];
					// Hibar3[k*N[0] + i] += Hitr3[k];
					// if(k <= 1){
					// 	Rprintf("i = %2i, Pi=%4.8e, Hibar2 = %4.8e, Hibar3 = %4.8e \n", i, Pi[i], Hibar2[k*N[0] + i], Hibar3[k*N[0] + i]);
					// }
				}
			} // end if_bhhh
			// end second term of the hessian
		} // end of loop within "r"
		// mean of Pir
		Pi[i] /= R[0];
		// Rprintf("i = %2i, Pi = %4.8f \n", i, Pi[i]);
		// calculate likelihood
		*lnls += log(Pi[i]);
		// calculate gradient
		for(k = 0; k < *K4; k++){
			grad[k] += gibar[k*N[0] + i];
		}
		// last term of the hessian: g_bar*g_bar' (no if *BHHH == 0 condition here)
		ijmiss = 0.0;
		for(jj = 0; jj < *K4; jj++){
			for(ii = 0; ii < *K4; ii++){
				if(ii >= jj){
					HibarBHHH[ jj * K4[0] + ii - ijmiss] = gibar[ii*N[0] + i] * gibar[jj*N[0] + i];
					// if(i <= 1){
					// 	Rprintf("i = %2i, HibarBHHH = %4.8E \n", i, HibarBHHH[ jj * K3[0] + ii - ijmiss]);
					// }
				} else {
					ijmiss++;
				}
			} //end "ii"
		}//end "jj"
		// calculate hessian
		if(*BHHH == 0){
			for(k = 0; k < k4H; k++){
				hesstri[k] += Hibar[k*N[0] + i] + Hibar2[k*N[0] + i] - HibarBHHH[k];
				// if(k == k4H-1){
				// 	Rprintf("i = %2i, H1 = %4.8E, H2 = %4.8E, H3 = %4.8E \n", i, Hibar[k*N[0] + i], Hibar2[k*N[0] + i], HibarBHHH[k]);
				// }
			}
		} // end if_bhhh
		// if BHHH
		if(*BHHH == 1){
			for(k = 0; k < k4H; k++){
				hesstri[k] -= HibarBHHH[k];
			}
		} // end if_bhhh
	} // end of loop within "i"
	
	// Rprintf("H_last = %4.8E \n", hesstri[k4H-1]);
		
	free(CZg);
	free(Pir);
	free(Pi);
	free(gitr);
	free(gir);
	free(gibar);
	free(Hitr);
	free(Hitr2);
	free(Hitr3);
	free(Hir);
	free(Hir2);
	// free(Hir3);
	free(Hibar);
	free(Hibar2);
	// free(Hibar3);
	free(HibarBHHH);
}

void gtre_ll(int *prod, double *W, double *H, int *N, int *R, double *Z, int *NT, double *C, 
double *ids, double *idvar, int *idlenmax, double *theta, int *K4,
double *lnls){
	
	// double lmd1 = exp( theta[*K3-3] );
	// double sig1 = exp( theta[*K3-2] );
	// double sigw = exp( theta[*K3-1] );
	// calculate CZg = C - Z*g
	// Allocate memory
	double * CZg   = (double *) malloc (*NT * sizeof(double));
	// Multiplication of Z and theta
	for(int q = 0; q < *NT; q++){
		for(int w = 0; w < *K4-4; w++){
			if(w == 0){
				CZg[q] = C[q];
			}
			CZg[q] -= Z[ w*NT[0] + q ] * theta[w];
		}
		// Rprintf("q=%2i, CZg=%4.3E \n", q, CZg[q]);
	}
	// Allocate memory
	double * Pir    = (double *) malloc (*R  * sizeof(double));
	double * Pi     = (double *) malloc (*N  * sizeof(double));	

	int r, i, j;
	// initialise those that come from R
	*lnls = 0.0;
	
	// get the sums of Pir first to get Qir later
	for(i = 0; i < *N; i++){
		// now within r
		// initialize for likelihood
		Pi[i] = 0.0;
		// go on
		for(r = 0; r < *R; r++){
			// now within i
			// int Ti = 0;
			// initialize for likelihood
			Pir[r] = 1.0;
			// go on
			for(j = 0; j < *NT; j++)
			{
				if(idvar[j] == ids[i])
				{
					// here is the main action
					double A = (CZg[j] - theta[*K4-2] * W[r*N[0] + i ] + theta[*K4-1] * H[r*N[0] + i ]*prod[0] ) / theta[*K4-3];
					double B = A * theta[*K4-4];
					// begin calculate "itr" element
					double Pitr = 2.0 * pow(theta[*K4-3],-1.0) * dnorm(A, 0, 1, 0) * pnorm(-B*prod[0], 0, 1, 1, 0);
					// end calculate "itr" element
					// product over "t"
					// Pir[r] *= Pitr;
					Pir[r] *= Pitr;
				} // end condition "idvar[j] == ids[i]"
			} // end checking if in this panel
			// sum of Pir
			Pi[i] += Pir[r];
		} // end of loop within "r"
		// Rprintf("i=%2i, Pi=%4.3E \n", i, Pi[i]);
		// mean of Pir
		Pi[i] /= R[0];
		// Rprintf("i = %2i, Pi = %4.8f \n", i, Pi[i]);
		// calculate likelihood
		*lnls += log(Pi[i]);
	} // end of loop within "i"
		
	free(CZg);
	free(Pir);
	free(Pi);
}


void gtre_grad(int *prod, double *W, double *H, int *N, int *R, double *Z, int *NT, double *C, 
double *ids, double *idvar, int *idlenmax, double *theta, int *K4,
double *lnls, double *grad){
	// int k4H = ( pow(K4[0],2) + K4[0])/2;
	// calculate CZg = C - Z*g
	// Allocate memory
	double * CZg   = (double *) malloc (*NT * sizeof(double));
	// Multiplication of Z and theta
	for(int q = 0; q < *NT; q++){
		for(int w = 0; w < *K4-4; w++){
			if(w == 0){
				CZg[q] = C[q];
			}
			CZg[q] -= Z[ w*NT[0] + q ] * theta[w];
		}
	}
	// Allocate memory
	double * Pir    = (double *) malloc (*R  * sizeof(double));
	double * Pi     = (double *) malloc (*N  * sizeof(double));
	double * gitr   = (double *) malloc (*K4 * sizeof(double));
	double * gir    = (double *) malloc (R[0]*K4[0] * sizeof(double));
	double * gibar  = (double *) malloc (N[0]*K4[0] * sizeof(double));
	
	int k, r, i, j;
	// get the sums of Pir first to get Qir later
	for(i = 0; i < *N; i++){
		// now within r
		// initialize for likelihood
		Pi[i] = 0.0;
		// go on
		for(r = 0; r < *R; r++){
			// now within i
			// int Ti = 0;
			// initialize for likelihood
			Pir[r] = 1.0;
			// go on
			for(j = 0; j < *NT; j++)
			{
				if(idvar[j] == ids[i])
				{
					// here is the main action
					double A = (CZg[j] - theta[*K4-2] * W[r*N[0] + i ] + theta[*K4-1] * H[r*N[0] + i ]*prod[0] ) / theta[*K4-3];
					double B = A * theta[*K4-4];
					// begin calculate "itr" element
					double Pitr = 2.0 * pow(theta[*K4-3],-1.0) * dnorm(A, 0, 1, 0) * pnorm(-B*prod[0], 0, 1, 1, 0);
					// end calculate "itr" element
					// product over "t"
					Pir[r] *= Pitr;
				} // end condition "idvar[j] == ids[i]"
			} // end checking if in this panel
			// sum of Pir
			Pi[i] += Pir[r];
		} // end of loop within "r"
		// Rprintf("i=%2i, Pi=%4.3E \n", i, Pi[i]);
	} // end of loop within "i"
	
	// Now the code for calculation of lnls, gr and hess
	
	// initialise those that come from R
	*lnls = 0.0;
	
	for(k = 0; k < *K4; k++){
		grad[k] = 0.0;
	}
	for(i = 0; i < *N; i++){
		// now within r
		// initialize for gradient
		for(k = 0; k < *K4; k++){
			gibar[k*N[0] + i] = 0.0;
		}
		// go on
		for(r = 0; r < *R; r++){
			// now within i
			// int Ti = 0;
			// initialize for likelihood
			Pir[r] = 1.0;
			// initialize for gradient
			for(k = 0; k < *K4; k++){
				gir[k*R[0] + r] = 0.0;	
			}
			// go on
			for(j = 0; j < *NT; j++)
			{
				if(idvar[j] == ids[i])
				{
					// here is the main action
					double A = (CZg[j] - theta[*K4-2] * W[r*N[0] + i ] + theta[*K4-1] * H[r*N[0] + i ]*prod[0] ) / theta[*K4-3];
					double B = A * theta[*K4-4];
					double C = exp( dnorm(-B, 0, 1, 1) - pnorm(-B*prod[0], 0, 1, 1, 1));
					// double D = B + C;
					double E = pow(theta[*K4-3],-1.0) * (A + theta[*K4-4]*C*prod[0]);
					// double F = -1.0 * pow(theta[*K3-2],-2.0) * (1 + pow(theta[*K3-3], 2.0) * C * D);
					// double G = -1.0 * C * (1.0 - B * D);
					// @@@ likelihood function
					// begin calculate "itr" element
					double Pitr = 2.0 * pow(theta[*K4-3],-1.0) * dnorm(A, 0, 1, 0) * pnorm(-B*prod[0], 0, 1, 1, 0);
					// end calculate "itr" element
					// Ti++;
					// if(i == 1){
					// 	Rprintf("r = %2i, j = %2i, i = %2i, Pitr = %4.8f \n", r, j, i, Pi);
					// }
					// product over "t"
					Pir[r] *= Pitr;
					// @@@ gradient
					// begin calculate "itr" element
					// w.r.t. gammas
					for(k = 0; k < *K4-4; k++){
						gitr[k] = E * Z[ k*NT[0] + j ];
					}
					// w.r.t. lambda
					gitr[*K4-4] = -C *prod[0] * A;
					// w.r.t. sigma
					gitr[*K4-3] = -1.0 * pow(theta[*K4-3],-1.0) + E * A;
					// w.r.t. sigma_w
					gitr[*K4-2] = E * W[r*N[0] + i ];
					// w.r.t. sigma_h
					gitr[*K4-1] = E * -H[r*N[0] + i ]*prod[0];
					// end calculate "itr" element
					// sum over "t"
					for(k = 0; k < *K4; k++){
						gir[k*R[0] + r] += gitr[k];	
					}
				} // end condition "idvar[j] == ids[i]"
			} // end checking if in this panel
			// weighted sums for gradient, sums of columns of (supposedly matrix) gir
			for(k = 0; k < *K4; k++){
				gibar[k*N[0] + i] += gir[k*R[0] + r] * Pir[r] / Pi[i];
			}
			// end second term of the hessian
		} // end of loop within "r"
		// mean of Pir
		Pi[i] /= R[0];
		// Rprintf("i = %2i, Pi = %4.8f \n", i, Pi[i]);
		// calculate likelihood
		*lnls += log(Pi[i]);
		// calculate gradient
		for(k = 0; k < *K4; k++){
			grad[k] += gibar[k*N[0] + i];
		}
	} // end of loop within "i"
		
	free(CZg);
	free(Pir);
	free(Pi);
	free(gitr);
	free(gir);
	free(gibar);
}
