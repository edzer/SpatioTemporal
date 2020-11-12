#include "R.h"
#include "Rinternals.h"
#include <R_ext/Lapack.h>

//calculate dot product of two vectors
//input vectors and vector length
SEXP C_dotProd(SEXP v1R, SEXP v2R){

  //extract input from R
  double *v1 = REAL(v1R);
  double *v2 = REAL(v2R);
  int n1 = length(v1R);
  int n2 = length(v2R);
  //loop variables
  int n, j; 
  //variables holding return matrix
  double *result;
  SEXP resultR;

  PROTECT(resultR = allocMatrix(REALSXP, 1, 1));
  result = REAL(resultR);

  //set n to smallest of n1,n2
  if(n1 < n2){
    n = n1;
  }else{
    n = n2;
  }
  result[0] = 0;
  for(j=0; j<n; j++){
    result[0] += v1[j]*v2[j];
  }
  
  UNPROTECT(1); /* resultR */
  return resultR;
}


//calculate sum of the squared elements of a vector (|x|^2)
//input vector and vector size
SEXP C_norm2(SEXP matR){

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
  for(j=0; j<n; j++){
    result[0] += mat[j]*mat[j];
  }
  
  UNPROTECT(1); /* resultR */
  return resultR;
}
