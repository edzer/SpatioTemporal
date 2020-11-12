#include "R.h"
#include "Rinternals.h"

/* COMMON HELPER FUNCITIONS */
int max_loc_F(int* vec, int n){
  int i, max_size=0;
  for(i=0; i<n; ++i){
    if( max_size < vec[i] ){
      max_size = vec[i];
    }
  }
  return max_size;
}

/* GENERAL ERROR CHECKING COMMON TO ALL BLOCK MATRIX FUNCTIONS */
void errorCheckingF(int length_F, int length_ind, int* loc_ind, int n_loc){
  if( length_F!=length_ind ){
    error("dim(F)[1] != length(loc.ind)");
  }
  if( max_loc_F(loc_ind, length_ind)>n_loc ){
    error("max(loc.ind) > number of locations");
  }
}

//calculates F'*X, inputs:
//  A vector/matrix X
//  The matrix F (for a given column each row corresponds to the value of the
//    temporal trend at the time of that observations, trends vary between columns)
//  number of different locations
//  Vector of locations, giving the location of each observation. Needed to know
//    which element of sigma.res that we should multiply elements in F with.
//returns a m*n_loc-by-dim(X)[2] matrix
SEXP C_calc_tFX(SEXP XR, SEXP FR, SEXP n_locR, SEXP locIndR){
  //extract input from R
  double *X = REAL(XR);
  double *F = REAL(FR);
  SEXP dimF = getAttrib(FR, R_DimSymbol);
  int n_obs = INTEGER(dimF)[0];
  int m = INTEGER(dimF)[1];
  int n_loc = *INTEGER(n_locR);
  int n_x = length(XR)/n_obs;
  int *locInd = INTEGER(locIndR);
  //loop variables
  int n_tot, j, i, i_m;
  //variables holding return matrix
  double *result;
  SEXP resultR;

  //error checking
  if( length(XR) % n_obs != 0){
    error("length(X) not multiple of dim(F)[1]");
  }
  errorCheckingF(n_obs, length(locIndR), locInd, n_loc);
    
  n_tot = n_loc*m;

  PROTECT(resultR = allocMatrix(REALSXP, n_tot, n_x));
  result = REAL(resultR);

  //INITIALIZE TO 0
  memset(result, 0, n_tot*n_x*sizeof(double));

  for(j=0; j<n_x; ++j)
    for(i_m=0; i_m<m; ++i_m)
      for(i=0; i<n_obs; ++i)
	result[locInd[i]-1 + i_m*n_loc + j*n_tot] += 
	  F[i+i_m*n_obs] * X[i+j*n_obs];
  
  UNPROTECT(1); //resultR
  return resultR;
}


//calculates F_i*X_i, inputs:
//  A matrix X_i
//  The vector F (one "column" of F, corresponds to the value of the
//    temporal trend at the time of that observations)
//  Vector of locations, giving the location of each observation. Needed to
//    know which element of X that we should multiply elements in F with.
//returns a column vector with (n_obs*m) elements
SEXP C_calc_F_part_X(SEXP XR, SEXP FR, SEXP locIndR){
  //extract input from R
  double *X = REAL(XR);
  SEXP dimX = getAttrib(XR, R_DimSymbol);
  int n_loc = INTEGER(dimX)[0];
  int pi = INTEGER(dimX)[1];
  double *F = REAL(FR);
  int n_obs = length(FR);
  int *locInd = INTEGER(locIndR);
  //loop variables
  int i, i_pi;
  //variables holding return matrix
  double *result;
  SEXP resultR;

  //error checking
  errorCheckingF(n_obs, length(locIndR), locInd, n_loc);

  PROTECT(resultR = allocMatrix(REALSXP, n_obs, pi));
  result = REAL(resultR);

  //INITIALIZE TO 0
  memset(result, 0, n_obs*pi*sizeof(double));

  for(i_pi=0; i_pi<pi; ++i_pi)
    for(i=0; i<n_obs; ++i)
      result[i + i_pi*n_obs] += F[i] * X[locInd[i]-1 + n_loc*i_pi];
  UNPROTECT(1); //resultR
  return resultR;
}


//calculates F'*inv(sigma.res)*F, inputs:
//  The matrix inv(sigma.res), ordinarily calculated by calls to
//    make_sigma_res, make_chol_block, and inv_chol_block
//  The matrix F (for a given column each row corresponds to the value of the
//    temporal trend at the time of that observations, trends vary between columns)
//  number of different locations
//  Vector of locations, giving the location of each observation. Needed to know
//    which element of sigma.res that we should multiply elements in F with.
//  Vector with the size of each block
//returned matrix will be square with side n_loc*m (same as sigmaB)
//  the returned matrix is band diagonal with bandwidth n_loc.
SEXP C_calc_tFXF(SEXP iSigmaR, SEXP FR, SEXP n_locR, SEXP locIndR,
		 SEXP block_sizesR){

  double *iSigma = REAL(iSigmaR);
  SEXP dimSigma = getAttrib(iSigmaR, R_DimSymbol);
  int n_obs = INTEGER(dimSigma)[0];
  double *F = REAL(FR);
  SEXP dimF = getAttrib(FR, R_DimSymbol);
  int m = INTEGER(dimF)[1];
  int n_loc = *INTEGER(n_locR);
  int *locInd = INTEGER(locIndR);
  int *block_sizes = INTEGER(block_sizesR);
  int n_blocks = length(block_sizesR);
  //loop variables
  int i, k, l, i_m, loc1, loc2, n_tot, tot_block_sizes=0;
  //return variables
  double *result, *iS_F;
  SEXP resultR;

  //error checking
  if( n_obs!=INTEGER(dimSigma)[1] ){
    error("'mat' not a square matrix");
  }
  if( n_obs!=INTEGER(dimF)[0] ){
    error("dim(F)[1] != dim(mat)[1]");
  }
  errorCheckingF(n_obs, length(locIndR), locInd, n_loc);
  //check total size of all blocks
  for(i=0; i<n_blocks; ++i){
    tot_block_sizes += block_sizes[i];
  }
  if( n_obs != tot_block_sizes ){
    error("sum(block.sizes) != dim(mat)[1]");
  }

  n_tot = n_loc*m;
  
  PROTECT(resultR = allocMatrix(REALSXP, n_tot, n_tot));
  result = REAL(resultR);

  //INITIALIZE TO 0
  memset(result, 0, n_tot*n_tot*sizeof(double));

  //temporary storage for inv(Sigma)*F
  iS_F = Calloc( n_obs*n_tot, double);
  memset(iS_F, 0, n_obs*n_tot * sizeof(double));

  //calculate inv(Sigma)*F
  for(i_m=0; i_m<m; ++i_m){
    loc1=0;
    for(i=0; i<n_blocks; ++i){
      for(k=0; k<block_sizes[i]; ++k){
	for(l=0; l<block_sizes[i]; ++l){
	  loc2 = locInd[loc1+l]-1;
	  iS_F[loc1+k + loc2*n_obs + i_m*n_obs*n_loc] +=
	    F[loc1+l+i_m*n_obs]*iSigma[loc1+k + (loc1+l)*n_obs];
	}
      }
      loc1 += block_sizes[i];
    }
  }
  //calculate F' * (inv(Sigma)*F)
  for(l=0; l<n_tot; ++l)
    for(k=0; k<m; ++k)
      for(i=0; i<n_obs; ++i){
	loc1 = locInd[i]-1;
	result[loc1+k*n_loc + l*n_tot] += F[i + k*n_obs] * iS_F[i + l*n_obs];
      }
  
  Free(iS_F);

  UNPROTECT(1); //resultR
  return resultR;
}
