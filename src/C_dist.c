#include "R.h"
#include "Rinternals.h"

/* FUNCTION THAT INTERFACE WITH R */
//Compute Euclidian distances between two sets of points
//inputs are:
//  matrix with first set of points
//  matrix with second set of points
//  indicator of symmetric matrix.
SEXP C_dist(SEXP coord1R, SEXP coord2R, SEXP symmetricR){

  //extract input from R
  double *coord1 = REAL(coord1R);
  double *coord2 = REAL(coord2R);
  int symmetric = *INTEGER(symmetricR);
  //sizes of the matrices
  SEXP c1Dim = getAttrib(coord1R, R_DimSymbol);
  SEXP c2Dim = getAttrib(coord2R, R_DimSymbol);
  int n1 = INTEGER(c1Dim)[0];
  int n2 = INTEGER(c2Dim)[0];
  int d = INTEGER(c1Dim)[1];
  //loop variables
  int i, j, k;
  double tmp;
  //variables holding return matrix
  double *result;
  SEXP resultR;

  //error checking
  if( d!=INTEGER(c2Dim)[1] ){
    error("coord1 and coord2 should have the same number of columns.");
  }
  if( symmetric && n1!=n2 ){
    error("symmetric assumes that coordinates are the same (length).");
  }
  
  //allocate data for return
  PROTECT(resultR = allocMatrix(REALSXP, n1, n2));
  result = REAL(resultR);
  memset(result, 0, n1*n2*sizeof(double));

  if( symmetric!=0 ){
    //compute sum_k( (x1_ik-x1_jk)^2 ), upper half only
    for(k=0; k<d; ++k){
      for(j=1; j<n1; ++j){
	for(i=0; i<j; ++i){
	  tmp = coord1[i+k*n1]-coord1[j+k*n1];
	  result[i+j*n1] += tmp*tmp;
	}
      }
    }
    //sqrt of all elements
    for(j=1; j<n1; ++j){
      for(i=0; i<j; ++i){
	result[i+j*n1] = sqrt(result[i+j*n1]);
      }
    }
    //fill in lower half of the matrix (diagonal=0)
    for(j=0; j<n1; ++j){
      for(i=j+1; i<n1; ++i){
	result[i+j*n1] = result[j+i*n1];
      }
    }
  }else{
    //compute sum_k( (x1_ik-x2_jk)^2 )
    for(k=0; k<d; ++k){
      for(j=0; j<n2; ++j){
	for(i=0; i<n1; ++i){
	  tmp = coord1[i+k*n1]-coord2[j+k*n2];
	  result[i+j*n1] += tmp*tmp;
	}
      }
    }
    //sqrt of all elements
    for(j=0; j<n2; ++j){
      for(i=0; i<n1; ++i){
	result[i+j*n1] = sqrt(result[i+j*n1]);
      }
    }
  }

  UNPROTECT(1); /* resultR */
  return resultR;
}
