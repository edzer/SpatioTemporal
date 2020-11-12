#include "R.h"
#include "Rinternals.h"
#include <R_ext/Lapack.h>
#include <string.h>

//declaration of external FORTRAN functions
extern void F77_CALL(dtrsl)(double*, int*, int*, double*, int*, int*);
extern void F77_CALL(dpodi)(double*, int*, int*, double*, int*);
extern void F77_CALL(dpofa)(double*, int*, int*, int*);


/* GENERAL ERROR CHECKING COMMON TO ALL BLOCK MATRIX FUNCTIONS */
void errorCheckingBlocks(SEXP matDim, int *block_sizes, int n_blocks){
  int i, tot_block_sizes = 0;
  //is matrix square?
  if( INTEGER(matDim)[0]!=INTEGER(matDim)[1] ){
    error("'mat' not a square matrix.");
  }
  //is sum(blocks) != dim(matrix)
  for(i=0; i<n_blocks; ++i){
    tot_block_sizes += block_sizes[i];
  }
  if( tot_block_sizes != INTEGER(matDim)[0] ){
    error("sum(block.sizes) != dim(mat)[1]");
  }
  return;
}

void errorCheckingBmat(int n, int length_x){
  if( length_x % n != 0 ){
    error("length(X) not multiple of dim(mat)[1]");
  }
}

/* COMMON HELPER FUNCITIONS */
int max_loc_block(int* vec, int n){
  int i, max_size=0;
  for(i=0; i<n; ++i){
    if( max_size < vec[i] ){
      max_size = vec[i];
    }
  }
  return max_size;
}

/* FUNCTION THAT INTERFACE WITH R */
//cholesky inverse of a block diagonal matrix using the diagonal structure.
//inputs are:
//  vector with the size of each block
//  maximum size of any block
SEXP C_make_chol_block(SEXP block_sizesR, SEXP matR){

  //extract input from R
  int *block_sizes = INTEGER(block_sizesR);
  int n_blocks = length(block_sizesR);
  double *mat = REAL(matR);
  SEXP matDim = getAttrib(matR, R_DimSymbol);
  int n_tot = INTEGER(matDim)[0];
  int max_size = 0;
  //loop variables
  int i, j, k, n_this, n_cum, info;
  int pos_def=1;
  //variables holding return matrix
  double *result, *block;
  SEXP resultR;

  //error checking
  errorCheckingBlocks(matDim, block_sizes, n_blocks);

  PROTECT(resultR = allocMatrix(REALSXP, n_tot, n_tot));
  result = REAL(resultR);

  if(n_blocks==1){ //don't need to copy matrix if only one block.
    memcpy(result, mat, n_tot*n_tot*sizeof(double));
    
    F77_CALL(dpofa)(result, &n_tot, &n_tot, &info);
    //info info!=0 then the matrix is not positive definite
    pos_def = (info==0);
    //set lower block to zero
    for (i=0; i < n_tot; ++i){
      for(j=i+1; j < n_tot; ++j){
	result[j + i*n_tot] = 0;
      }
    }
  }else{
    //max block size
    max_size = max_loc_block(block_sizes, n_blocks);
    block = Calloc(max_size*max_size, double);
    memset(result, 0, n_tot*n_tot*sizeof(double));
    
    n_cum = 0;
    for(i=0; i < n_blocks; i++){
      n_this = block_sizes[i];
      for(k=0; k < n_this; ++k){
	for(j=0; j <=k; ++j){
	  block[j+n_this*k] = mat[j+n_cum + (k+n_cum)*n_tot];
	}
      }
      F77_CALL(dpofa)(block, &n_this, &n_this, &info);
      //info info!=0 then the matrix is not positive definite
      pos_def = (info==0);
      for (k=0; k < n_this; ++k){
	for(j=0; j <= k; ++j){
	  result[j+n_cum + (k+n_cum)*n_tot] = block[j+n_this*k]; 
	}
      }    
      n_cum += n_this;
    }
    Free(block);
  }
  if(!pos_def) //not positive definite, set top left value to -1
    result[0 + 0*n_tot] = -1;

  UNPROTECT(1); /* resultR */
  return resultR;
}


//calculate a mtrix inverse from a cholesky factor (requires a previous call
//to make_chol_block)
//inputs are:
//  vector with the size of each block
//  the matrix it self.
SEXP C_inv_chol_block(SEXP block_sizesR, SEXP matR){

  //extract input from R
  int *block_sizes = INTEGER(block_sizesR);
  int n_blocks = length(block_sizesR);
  double *mat = REAL(matR);
  SEXP matDim = getAttrib(matR, R_DimSymbol);
  int n_tot = INTEGER(matDim)[0];
  int max_size = 0;
  //loop variables
  int i, j, k, n_this, n_cum;
  int job=1; //compute inverse only (no determinant)
  //variables holding return matrix
  double *result, *block, det;
  SEXP resultR;

  //error checking
  errorCheckingBlocks(matDim, block_sizes, n_blocks);
  
  PROTECT(resultR = allocMatrix(REALSXP, n_tot, n_tot));
  result = REAL(resultR);

  if(n_blocks==1){ //don't need to copy matrix if only one block.
    memcpy(result, mat, n_tot*n_tot*sizeof(double));
	
    F77_CALL(dpodi)(result, &n_tot, &n_tot, &det, &job);
    //fill in symmetry
    for (i=1; i < n_tot; ++i){
      for(j=0; j < i; ++j){
	result[i + j*n_tot] = result[j + i*n_tot];
      }
    }
  }else{
    //max block size
    max_size = max_loc_block(block_sizes, n_blocks);
    block = Calloc(max_size*max_size, double);
    memset(result, 0, n_tot*n_tot*sizeof(double));
	
    n_cum = 0;
    for(i=0; i < n_blocks; ++i){
      n_this = block_sizes[i];
      for (k=0; k < n_this; ++k){
	for(j=0; j <=k; ++j){
	  block[j+n_this*k] = mat[j+n_cum + (k+n_cum)*n_tot];
	}
      }
      F77_CALL(dpodi)(block, &n_this, &n_this, &det, &job);
      for (k=0; k < n_this; ++k){
	for(j=0; j <= k; ++j){
	  result[j+n_cum + (k+n_cum)*n_tot] = block[j+n_this*k]; 
	  result[k+n_cum + (j+n_cum)*n_tot] = block[j+n_this*k]; 
	}
      }
      n_cum += n_this;
    }      
    Free(block);
  }
  
  UNPROTECT(1); /* resultR */
  return resultR;
}

//solves the system T*X=B where T is block upper triagular (output from
//  make_chol_block)
//inputs is:
//  vector with the size of each block
//  indicator wether or not to transpose T
//  the triangular matrix
//  the LHS (B) of the equation system.
SEXP C_solve_tri_block(SEXP block_sizesR, SEXP transposeR, SEXP matR, SEXP xR){

  //extract input from R
  int *block_sizes = INTEGER(block_sizesR);
  int n_blocks = length(block_sizesR);
  double *mat = REAL(matR);
  SEXP matDim = getAttrib(matR, R_DimSymbol);
  int n_tot = INTEGER(matDim)[0];
  double *x = REAL(xR);
  int n_x = length(xR)/n_tot;
  int max_size = 0;
  //loop variables
  int i, j, k, n_this, n_cum,info, job;
  //variables holding return matrix
  double *result, *block;
  SEXP resultR;

  //error checking
  errorCheckingBlocks(matDim, block_sizes, n_blocks);
  errorCheckingBmat(n_tot, length(xR));
  //max block size
  max_size = max_loc_block(block_sizes, n_blocks);
  
  //job defines T*X=B to solve. T is upper triangular (=1)
  //and possibly transpose (+10)
  job = 10*(*INTEGER(transposeR)!=0)+1; 
  PROTECT(resultR = allocMatrix(REALSXP, n_tot, n_x));
  result = REAL(resultR);

  memcpy(result, x, n_tot*n_x*sizeof(double));

  if(n_blocks==1){ //don't need to copy matrix if only one block.
    for(i=0; i<n_x; ++i){
      F77_CALL(dtrsl)(mat, &n_tot, &n_tot, result+i*n_tot, &job, &info);
    }
  }else{
    block = Calloc(max_size*max_size, double);
    n_cum = 0;
    for(k=0; k < n_blocks; ++k){
      n_this = block_sizes[k];
      for (i=0; i < n_this; ++i){
	for(j=0; j <=i; ++j){
	  block[j+n_this*i] = mat[j+n_cum + (i+n_cum)*n_tot];
	}
      }
      for(i=0; i<n_x; ++i){
	F77_CALL(dtrsl)(block, &n_this, &n_this, result+n_cum+i*n_tot, &job, &info);
      }
      n_cum += n_this;
    }
    Free(block);
  }
  UNPROTECT(1); /* resultR */
  return resultR;
}


//multiplication of block diagonal matrix with another matrix
//inputs are:
//  vector with the size of each block
//  the block-matrix it self
//  the matrix to multiply with
// BLOCK_MATRIX*X
SEXP C_block_mult(SEXP block_sizesR, SEXP matR, SEXP xR){

  //extract input from R
  int *block_sizes = INTEGER(block_sizesR);
  int n_blocks = length(block_sizesR);
  double *mat = REAL(matR);
  SEXP matDim = getAttrib(matR, R_DimSymbol);
  int n_tot = INTEGER(matDim)[0];
  double *x = REAL(xR);
  int n_x = length(xR)/n_tot;
  //loop variables
  int i, k, l, j, n_cum;
  //variables holding return matrix
  double *result;
  SEXP resultR;

  //error checking
  errorCheckingBlocks(matDim, block_sizes, n_blocks);
  errorCheckingBmat(n_tot, length(xR));
  
  PROTECT(resultR = allocMatrix(REALSXP, n_tot, n_x));
  result = REAL(resultR);

  //set results to zero
  memset(result, 0, n_tot*n_x*sizeof(double));
    
  for(j=0; j<n_x; ++j){
    n_cum = 0;
    for(i=0; i < n_blocks; ++i){
      for(k=0; k < block_sizes[i]; ++k)
	for(l=0; l < block_sizes[i]; ++l)
	  result[n_cum+k+j*n_tot] += mat[n_cum+k + (n_cum+l)*n_tot] * 
	    x[n_cum+l+j*n_tot];
      n_cum += block_sizes[i];
    }
  }

  UNPROTECT(1); /* resultR */
  return resultR;
}
