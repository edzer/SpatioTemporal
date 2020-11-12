#include "R.h"
#include "Rinternals.h"
#include <R_ext/Lapack.h>

//calculate sum of the logarithm of diagonal elements, corresponds to the
//log determinant of a cholesky factor
//input matrix
SEXP C_sum_log_diag(SEXP matR){

  //extract input from R
  double *mat = REAL(matR);
  SEXP matDim = getAttrib(matR, R_DimSymbol);
  int n = INTEGER(matDim)[0];
  //loop variables
  int j;
  //variables holding return matrix
  double *result;
  SEXP resultR;

  if( n!=INTEGER(matDim)[1] ){
    error("'mat' not a square matrix.");
  }
  PROTECT(resultR = allocMatrix(REALSXP, 1, 1));
  result = REAL(resultR);

  result[0] = 0;
  for(j=0; j<n; j++)
    result[0] += log(mat[j+n*j]);
  
  UNPROTECT(1); /* resultR */
  return resultR;
}


//calculate sum of the log elements of a vector; sum(log(x))
//input vector and vector size
SEXP C_sum_log(SEXP matR){

  //extract input from R
  double *mat = REAL(matR);
  int n = length(matR);
  //loop variables
  int j;
  //variables holding return matrix
  double *result;
  SEXP resultR;

  PROTECT(resultR = allocMatrix(REALSXP, 1, 1));
  result = REAL(resultR);

  result[0] = 0;
  for(j=0; j<n; j++)
    result[0] += log(mat[j]);
  
  UNPROTECT(1); /* resultR */
  return resultR;
}
