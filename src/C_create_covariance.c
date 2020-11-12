#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"
#include <string.h>


/* ADDING A NEW COVARIANCE FUNCTION */
//1) create a covariance function, see exp_cov etc below. IMPORTANT NUGGET IS
//   ADDED AT A LATER STAGE!
//2) create a function that returns the names of parameters, existing prototypes are:
//     SEXP cov_names_none()
//     SEXP cov_names_range_sill()
//     SEXP cov_names_range_sill_shape()
//3) Register your functions (with suitable name) in
//     cov_funs_struct selectCov(const char* cov_type)
//4) Add the name of your covariance function to
//     SEXP C_cov_names()
//5) Add documentation to the R function
//     namesCovFuns
//   in
//     c_cov_matricies.R


/* EXPONENTIAL COVARIANCE FUNCTION */
//given inputs n, and pointers to parameters and distance matrix
//it will compute the covariance function for all n elements in dist and
//return them in the PREALLOCATED array cov.
void exp_cov(int n, double* par, double* dist, double* cov, int* diff){
  //extract parameters
  int i;
  double range = par[0];
  double sill = par[1];
  double C1 = 0;
  if( diff[0]==0 && diff[1]==0 ){
    //f(d)
    for(i=0; i<n; i++){
      cov[i] = sill*exp(-dist[i]/range);
    }
  }else if(diff[1]==0){
    if( diff[0]==1 ){
      //f'_range(d)
      C1 = sill/R_pow_di(range, 2);
      for(i=0; i<n; i++){
        cov[i] = C1*dist[i]*exp(-dist[i]/range);
      }
    }else{
      //f'_sill(d)
      for(i=0; i<n; i++){
        cov[i] = exp(-dist[i]/range);
      }
    }
  }else{
    if( diff[0]==1 && diff[1]==1 ){
      //f''_range(d)
      C1 = sill/R_pow_di(range, 4);
      for(i=0; i<n; i++){
        cov[i] = exp(-dist[i]/range) * C1 * dist[i] * (dist[i]-2*range);
      }
    }else if( diff[0]==2 && diff[1]==1 ){
      //f''_range,sill(d)
      C1 = R_pow_di(range, 2);
      for(i=0; i<n; i++){
        cov[i] = dist[i]*exp(-dist[i]/range)/C1;
      }
    }else{
      //f''_sill(d)
      memset( cov, 0, n*sizeof(double) );
    }
  }
  return;
}

/* DOUBLE EXPONENTIAL/GAUSSIAN COVARIANCE FUNCTION */
void exp2_cov(int n, double* par, double* dist, double* cov, int* diff){
  //extract parameters
  int i;
  double range = par[0];
  double range2 = range*range;
  double sill = par[1];
  double C1=0;
  if( diff[0]==0 && diff[1]==0 ){
    //f(d)
    for(i=0; i<n; i++){
      cov[i] = sill*exp(-dist[i]*dist[i]/range2);
    }
  }else if(diff[1]==0){
    if( diff[0]==1 ){
      //f'_range(d)
      C1 = 2*sill/R_pow_di(range, 3);
      for(i=0; i<n; i++){
        cov[i] = C1*dist[i]*dist[i]*exp(-dist[i]*dist[i]/range2);
      }
    }else{
      //f'_sill(d)
      for(i=0; i<n; i++){
        cov[i] = exp(-dist[i]*dist[i]/range2);
      }
    }
  }else{
    if( diff[0]==1 && diff[1]==1 ){
      //f''_range(d)
      C1 = 2*sill/R_pow_di(range, 6);
      for(i=0; i<n; i++){
        cov[i] = exp(-dist[i]*dist[i]/range2) * C1*dist[i]*dist[i] *
          (2*dist[i]*dist[i] - 3*range2);
      }
    }else if( diff[0]==2 && diff[1]==1 ){
      //f''_range,sill(d)
      C1 = 2/(range2*range);
      for(i=0; i<n; i++){
        cov[i] = C1*dist[i]*dist[i]*exp(-dist[i]*dist[i]/range2);
      }
    }else{
      //f''_sill(d)
      memset( cov, 0, n*sizeof(double) );
    }
  }

  return;
}

/* CUBIC COVARIANCE FUNCTION */
void cubic_cov(int n, double* par, double* dist, double* cov, int* diff){
  //extract parameters
  int i;
  double range = par[0];
  double sill = par[1];
  double frac, frac2, frac3, frac5, frac7;
  memset( cov, 0, n*sizeof(double) );

  if( diff[0]==0 && diff[1]==0 ){
    //f(d)
    for(i=0; i<n; i++){
      frac = dist[i]/range;
      if( frac<1 ){
        frac2 = R_pow_di(frac, 2);
        frac3 = R_pow_di(frac, 3);
        frac5 = R_pow_di(frac, 5);
        frac7 = R_pow_di(frac, 7);
        cov[i] = 1.0 - (7.0*frac2 - 8.75*frac3 + 3.5*frac5 - 0.75*frac7);
        cov[i] *= sill;
      }
    }
  }else if(diff[1]==0){
    if( diff[0]==1 ){
      //f'_range(d)
      for(i=0; i<n; i++){
        frac = dist[i]/range;
        if( frac<1 ){
          frac2 = R_pow_di(frac, 2)/range;
          frac3 = R_pow_di(frac, 3)/range;
          frac5 = R_pow_di(frac, 5)/range;
          frac7 = R_pow_di(frac, 7)/range;
          cov[i] = 14.0*frac2 - 26.25*frac3 + 17.5*frac5 - 5.25*frac7;
          cov[i] *= sill;
        }
      }
    }else{
      //f'_sill(d)
      for(i=0; i<n; i++){
        frac = dist[i]/range;
        if( frac<1 ){
          frac2 = R_pow_di(frac, 2);
          frac3 = R_pow_di(frac, 3);
          frac5 = R_pow_di(frac, 5);
          frac7 = R_pow_di(frac, 7);
          cov[i] = 1.0 - (7.0*frac2 - 8.75*frac3 + 3.5*frac5 - 0.75*frac7);
        }
      }
    }
  }else{
    if( diff[0]==1 && diff[1]==1 ){
      //f''_range(d)
      for(i=0; i<n; i++){
        frac = dist[i]/range;
        if( frac<1 ){
          frac2 = R_pow_di(frac, 2) / R_pow_di(range, 2);
          frac3 = R_pow_di(frac, 3) / R_pow_di(range, 2);
          frac5 = R_pow_di(frac, 5) / R_pow_di(range, 2);
          frac7 = R_pow_di(frac, 7) / R_pow_di(range, 2);
          cov[i] = -42.0*frac2 + 105.0*frac3 - 105.0*frac5 + 42.0*frac7;
          cov[i] *= sill;
        }
      }
    }else if( diff[0]==2 && diff[1]==1 ){
      //f''_range,sill(d)
      for(i=0; i<n; i++){
        frac = dist[i]/range;
        if( frac<1 ){
          frac2 = R_pow_di(frac, 2)/range;
          frac3 = R_pow_di(frac, 3)/range;
          frac5 = R_pow_di(frac, 5)/range;
          frac7 = R_pow_di(frac, 7)/range;
          cov[i] = 14.0*frac2 - 26.25*frac3 + 17.5*frac5 - 5.25*frac7;
        }
      }
    }else{
      //f''_sill(d), zero second derivative already set above
    }
  }
  
  return;
}

/* SPHERICAL COVARIANCE FUNCTION */
void spherical_cov(int n, double* par, double* dist, double* cov, int* diff){
  //extract parameters
  int i;
  double range = par[0];
  double sill = par[1];
  double frac, frac3;
  memset( cov, 0, n*sizeof(double) );
  if( diff[0]==0 && diff[1]==0 ){
    //f(d)
    for(i=0; i<n; i++){
      frac = dist[i]/range;
      if( frac<1 ){
        frac3 = R_pow_di(frac, 3);
        cov[i] = 1.0 - 1.5*frac + 0.5*frac3;
        cov[i] *= sill;
      }
    }
  }else if(diff[1]==0){
    if( diff[0]==1 ){
      //f'_range(d)
      for(i=0; i<n; i++){
        frac = dist[i]/range;
        if( frac<1 ){
          frac3 = R_pow_di(frac, 3)/range;
          frac /= range;
          cov[i] = 1.5*frac - 1.5*frac3;
          cov[i] *= sill;
        }
      }
    }else{
      //f'_sill(d)
      for(i=0; i<n; i++){
        frac = dist[i]/range;
        if( frac<1 ){
          frac3 = R_pow_di(frac, 3);
          cov[i] = 1.0 - 1.5*frac + 0.5*frac3;
        }
      }
    }
  }else{
    if( diff[0]==1 && diff[1]==1 ){
      //f''_range(d)
      for(i=0; i<n; i++){
        frac = dist[i]/range;
        if( frac<1 ){
          frac3 = R_pow_di(frac, 3) / R_pow_di(range, 2);
          frac /= R_pow_di(range, 2);
          cov[i] = -3.0*frac + 6.0*frac3;
          cov[i] *= sill;
        }
      }
    }else if( diff[0]==2 && diff[1]==1 ){
      //f''_range,sill(d)
      for(i=0; i<n; i++){
        frac = dist[i]/range;
        if( frac<1 ){
          frac3 = R_pow_di(frac, 3)/range;
          frac /= range;
          cov[i] = 1.5*frac - 1.5*frac3;
        }
      }
    }else{
      //f''_sill(d), zero second derivative already set above
    }
  }

  return;
}

/* MATERN COVARIANCE FUNCTION */
void matern_cov(int n, double* par, double* dist, double* cov, int* diff){
  //extract parameters
  int i;
  double sill = par[1];
  double shape = par[2];
  //scaled range
  double scaled_range = sqrt(8*shape)/par[0];
  //log of normalising constant sill / (2^(nu-1)*gamma(nu))
  double log_norm = log(sill) - (shape-1)*log(2) - lgammafn(shape);
  double tmp, log_value;

  if( diff[0]==0 && diff[1]==0 ){
    //f(d)
    for(i=0; i<n; i++){
      if( dist[i]==0 ){
        cov[i] = sill;
      }else{
        tmp = scaled_range*dist[i];
        //recall that R bessel returns exp(x)*K_nu(x)
        //so we want log( exp(x)*K_nu(x)*exp(-x) ) = log( exp(x)*K_nu(x) ) - x
        log_value = log_norm + shape*log(tmp) + log(bessel_k(tmp, shape, 2)) - tmp;
        //transform back to non-log scale.
        cov[i] = exp( log_value );
      }
    }
  }else{
    error("Analytical derivatives not available for Matern covariance.");
  }
  return;
}

/* CAUCHY COVARIANCE FUNCTION */
void cauchy_cov(int n, double* par, double* dist, double* cov, int* diff){
  //extract parameters
  int i;
  double range = par[0];
  double range2 = range*range;
  double sill = par[1];
  double shape = par[2];
  double C1 = 0;
  double C2 = 0;

  if( diff[0]==0 && diff[1]==0 ){
    //f(d)
    for(i=0; i<n; i++){
      cov[i] = sill * R_pow(1 + dist[i]*dist[i]/range2, -shape);
    }
  }else if(diff[1]==0){
    if( diff[0]==1 ){
      //f'_range(d)
      C1 = 2*sill*shape / R_pow_di(range,3);
      for(i=0; i<n; i++){
        cov[i] =  C1 * dist[i]*dist[i] * R_pow(1 + dist[i]*dist[i]/range2, -shape-1);
      }
    }else if( diff[0]==2 ){
      //f'_sill(d)
      for(i=0; i<n; i++){
        cov[i] = R_pow(1 + dist[i]*dist[i]/range2, -shape);
      }
    }else{
      //f'_shape(d)
      for(i=0; i<n; i++){
        cov[i] = -sill * R_pow(1 + dist[i]*dist[i]/range2, -shape) *
          log(1 + dist[i]*dist[i]/range2);
      }
    }
  }else{
    if( diff[0]==1 && diff[1]==1 ){
      //f''_range(d)
      C1 = 4*sill*shape*(shape+1) / R_pow_di(range,6);
      C2 = -6*sill*shape / R_pow_di(range,4);
      for(i=0; i<n; i++){
        cov[i] =  C1 * R_pow_di(dist[i], 4) *
          R_pow(1 + dist[i]*dist[i]/range2, -shape-2) +
          C2 * dist[i]*dist[i] * R_pow(1 + dist[i]*dist[i]/range2, -shape-1);
      }
    }else if( diff[0]==2 && diff[1]==1 ){
      //f''_range,sill(d)
      C1 = 2*shape / R_pow_di(range,3);
      for(i=0; i<n; i++){
        cov[i] =  C1 * dist[i]*dist[i] * R_pow(1 + dist[i]*dist[i]/range2, -shape-1);
      }
    }else if( diff[0]==3 && diff[1]==1 ){
      //f''_range,shape(d)
      C1 = 2*sill / R_pow_di(range,3);
      for(i=0; i<n; i++){
        cov[i] =  C1 * dist[i]*dist[i] *
          R_pow(1 + dist[i]*dist[i]/range2, -shape-1) *
          (1 - shape * log(1 + dist[i]*dist[i]/range2));
      }
    }else if( diff[0]==2 && diff[1]==2 ){
      //f''_sill(d), zero second derivative zero
      memset( cov, 0, n*sizeof(double) );
    }else if( diff[0]==3 && diff[1]==2 ){
      //f''_sill,shape(d)
      for(i=0; i<n; i++){
        cov[i] = - R_pow(1 + dist[i]*dist[i]/range2, -shape) *
          log(1 + dist[i]*dist[i]/range2);
      }
    }else{
      //f''_shape,shape(d)
      for(i=0; i<n; i++){
        cov[i] = sill * R_pow(1 + dist[i]*dist[i]/range2, -shape) *
          log(1 + dist[i]*dist[i]/range2) * log(1 + dist[i]*dist[i]/range2);
      }
    }
  }
  
  return;
}

/* IID COVARIANCE FUNCTION (0 matrix, since nugget is added later*/
void iid_cov(int n, double* par, double* dist, double* cov, int* diff){
  memset( cov, 0, n*sizeof(double) );
  return;
}

/* COMMON SETS OF NAMES FOR COVARIANCE PARAMETERS */
//no parameters (e.g. iid)
SEXP cov_names_none(){
  SEXP parnames;           
  PROTECT( parnames=allocVector(STRSXP,0) );
  UNPROTECT(1);
  return parnames;
}

//just range and sill (e.g. exp, exp2/gaussian, cubic, spherical
SEXP cov_names_range_sill(){
  SEXP parnames;           
  PROTECT( parnames=allocVector(STRSXP,2) );
  SET_STRING_ELT(parnames, 0, mkChar("range"));
  SET_STRING_ELT(parnames, 1, mkChar("sill"));
  UNPROTECT(1);
  return parnames;
}

//range, sill, and shape (e.g. matern, cauchy)
SEXP cov_names_range_sill_shape(){
  SEXP parnames;           
  PROTECT( parnames=allocVector(STRSXP,3) );
  SET_STRING_ELT(parnames, 0, mkChar("range"));
  SET_STRING_ELT(parnames, 1, mkChar("sill"));
  SET_STRING_ELT(parnames, 2, mkChar("shape"));
  UNPROTECT(1);
  return parnames;
}

/* SELECTION OF COVARIANCE FUNCTION */
//First lets define a struct containing the two function pointers
typedef struct 
{
  void (*cov)(int, double*, double*, double*, int*);
  SEXP (*cov_names)();
} cov_funs_struct;

cov_funs_struct selectCov(const char* cov_type){
  //declare structure containing function pointers
  cov_funs_struct cov_funs = {NULL, NULL};
  //figure out which covariance function we want
  if( strcmp(cov_type, "exp")==0 || strcmp(cov_type, "exponential")==0 ){
    //exponential
    cov_funs.cov = exp_cov;
    cov_funs.cov_names = cov_names_range_sill;
  }else if( strcmp(cov_type, "exp2")==0 || strcmp(cov_type, "exponential2")==0 ||
            strcmp(cov_type, "gaussian")==0 ){
    //double exponential / Gaussian
    cov_funs.cov = exp2_cov;
    cov_funs.cov_names = cov_names_range_sill;
  }else if( strcmp(cov_type, "cubic")==0 ){
    //cubic
    cov_funs.cov = cubic_cov;
    cov_funs.cov_names = cov_names_range_sill;
  }else if( strcmp(cov_type, "spherical")==0 ){
    //spherical
    cov_funs.cov = spherical_cov;
    cov_funs.cov_names = cov_names_range_sill;
  }else if( strcmp(cov_type, "matern")==0 ){
    //matern
    cov_funs.cov = matern_cov;
    cov_funs.cov_names = cov_names_range_sill_shape;
  }else if( strcmp(cov_type, "cauchy")==0 ){
    //cauchy
    cov_funs.cov = cauchy_cov;
    cov_funs.cov_names = cov_names_range_sill_shape;
  }else if( strcmp(cov_type, "iid")==0 ){
    //iid
    cov_funs.cov = iid_cov;
    cov_funs.cov_names = cov_names_none;
  }
  return cov_funs;
}

/* (INLINE) HELPER FUNCTIONS FOR SIGMA_B AND SIGMA_NU */
//function that computes the covariance for one block in the block matrix
//inputs are:
//  function pointers
//  pointer to parameters
//  random effect variance to add to all elements
//  pointer to distance matrix
//  pointer to the matrix holding the covariance function.
//  size of distance matrix n1-x-n2
//  offset after each column in the distance matrix, (n.blocks-1)*n1
//  indicator if the matrix is symmetric.
void computeCovarianceBlock(cov_funs_struct cov_funs, double* pars,
                            double random, double* dist, double* covf,
                            int n1, int n2, int n_offset, int symmetric,
                            int* diff){
  //loop variables.
  int i, j;
  //total offset is n1+n_offset
  int n_offset_tot = n_offset+n1;
  
  //three different cases:
  if( symmetric!=0 && n1==n2 ){
    //1) symmetric (also require n1==n2)
    //compute covariance for the upper right half if the i:th column
    for(i=0; i<n2; ++i){
      cov_funs.cov(i+1, pars, dist+i*n1, covf+i*n_offset_tot, diff);
    }
    //Use symmetry to fill in the lower left half
    for(i=0; i<n2; ++i){
      for(j=i+1; j<n1; ++j){
        covf[j+i*n_offset_tot] = covf[i+j*n_offset_tot];
      }
    }
  }else if( n_offset!=0 ){
    //2) non-symmetric, n_offset!=0
    //compute covariance for the i:th column
    for(i=0; i<n2; ++i){
      cov_funs.cov(n1, pars, dist+i*n1, covf+i*n_offset_tot, diff);
    }
  }else{
    //3) non-symmetric, n_offset==0 
    //compute everything at once.
    cov_funs.cov(n1*n2, pars, dist, covf, diff);
  }//if( symmetric!=0 && n1==n2 ){...}else if( n_offset!=0 ){...}else{...}
  
  //if random effect is non-zero, add to all elements in the block
  if( random!=0 ){
    for(i=0; i<n2; ++i){
      for(j=0; j<n1; ++j){
        covf[j+i*n_offset_tot] += random;
      }
    }
  }//if( random!=0 )
  return;
}

/* R INTERFACE FUNCTIONS FOR TYPE OF COVARIANCE/PARAMETER NAMES*/
//first, return possible names of covariance functions
SEXP C_cov_names(){
  SEXP parnames;
  //increase the number to account for additional covariance names
  int n_covs = 10;
  PROTECT( parnames=allocVector(STRSXP,n_covs) );
  //add name(s) to the end of the list
  SET_STRING_ELT(parnames, 0, mkChar("exp"));
  SET_STRING_ELT(parnames, 1, mkChar("exponential"));
  SET_STRING_ELT(parnames, 2, mkChar("exp2"));
  SET_STRING_ELT(parnames, 3, mkChar("exponential2"));
  SET_STRING_ELT(parnames, 4, mkChar("gaussian"));
  SET_STRING_ELT(parnames, 5, mkChar("cubic"));
  SET_STRING_ELT(parnames, 6, mkChar("spherical"));
  SET_STRING_ELT(parnames, 7, mkChar("matern"));
  SET_STRING_ELT(parnames, 8, mkChar("cauchy"));
  SET_STRING_ELT(parnames, 9, mkChar("iid"));
   
  //release memory back to R
  UNPROTECT(1);
  return parnames;
}

//second, given a string, return expected names of parameters
SEXP C_cov_pars(SEXP cov_typeR){

  const char* cov_type;
  cov_funs_struct cov_funs;
  if( length(cov_typeR)!=1 ){
    error("'type' has to be length=1, is %d",length(cov_typeR));
  }

  //get name from R variable
  cov_type = CHAR( STRING_ELT(cov_typeR,0) );
  //figure out the matching function pointers
  cov_funs = selectCov(cov_type);

  if( cov_funs.cov==NULL || cov_funs.cov_names==NULL ){
    //bad name of covariance function, return NULL
    return R_NilValue;
  }
  //else return names assigned above
  return cov_funs.cov_names();
}

/* COMMON ERROR CHECKING FUNCTIONS */
int max_create_cov(int* vec, int n){
  int i, max_size=0;
  for(i=0; i<n; ++i){
    if( max_size < vec[i] ){
      max_size = vec[i];
    }
  }
  return max_size;
}

void checkCovarianceFunAndPars(cov_funs_struct cov_funs, const char* cov_type,
                               int n_pars, int* diff){
  //valid covariance function?
  if( cov_funs.cov==NULL || cov_funs.cov_names==NULL ){
    error("Unknown covariance specification: 'type' = %s",
          cov_type);
  }
  //check length of parameter vector
  if( length(cov_funs.cov_names())!=n_pars ){
    error("Expected %d parameter(s), but length(pars) = %d",
          n_pars, length(cov_funs.cov_names()));
  }
  if( diff!=NULL ){ //TODO, drop NULL check
    //check that diff parameters are reasonable
    if( diff[0]<0 || diff[0]>n_pars || diff[1]<0 || diff[1]>n_pars ){
      error("Elements in diff must be between 0 and %d; are: %d,%d",
            n_pars, diff[0], diff[1]);
    }
  }
}//checkCovarianceFunAndPars(...)

/* COMMON HELPER FUNCTIONS AND STRUCTS */
//extract information regarding differentiation (into pre allocated 2-element vector)
void extractDiffInfo(SEXP diffR, int* diff){
  diff[0] = INTEGER(diffR)[0];
  if( length(diffR)>1 ){
    int d1 = INTEGER(diffR)[1];
    //ensure that diff[0] >= diff[1]
    if( d1 > diff[0] ){
      diff[1] = diff[0];
      diff[0] = d1;
    }else{
      diff[1] = d1;
    }
  }else{
    diff[1] = 0;
  }
}


/* R INTERFACE FUNCTIONS FOR COMPUTING COVARIANCES */
//simple covariance function that just computes values for
//all elements in the given distance matrix/vector
SEXP C_cov_simple(SEXP cov_typeR, SEXP parsR, SEXP distR, SEXP diffR){

  //get covariance function name from R variable
  const char* cov_type = CHAR( STRING_ELT(cov_typeR,0) );
  //pointer to parameters and distance matrix
  double *pars = REAL(parsR);
  double *dist = REAL(distR);
  //number of elements in distance matrix
  int n = length(distR);
  int n_pars = length(parsR);
  //extract differentiation flags
  int diff[2];
  extractDiffInfo(diffR, diff);
  
  //return variable
  SEXP resultR;
  //figure out the function pointers
  cov_funs_struct cov_funs = selectCov(cov_type);

  //check if covariance function is valid
  checkCovarianceFunAndPars(cov_funs, cov_type, n_pars, diff);
  
  //allocate memory for return value
  PROTECT( resultR = allocMatrix(REALSXP, n, 1) );

  //compute covariance function
  cov_funs.cov(n, pars, dist, REAL(resultR), diff);

  //return values
  UNPROTECT(1); /* resultR */
  return resultR;
}//C_cov_simple(...)

//compute block cross-covariance matrix for Sigma_B
SEXP C_makeSigmaB(SEXP cov_typeR, SEXP parsR, SEXP nuggetR, SEXP distR,
                  SEXP symmetryR, SEXP ind2_to_1R, SEXP sparseR, SEXP diffR){
  //EXTRACT PARAMETERS FROM R CALL
  int n_blocks = length(cov_typeR);
  double* nugget = REAL(nuggetR);
  double* dist = REAL(distR);
  //dimensions of the distance matrix
  SEXP matDim = getAttrib(distR, R_DimSymbol);
  int n1 = INTEGER(matDim)[0];
  int n2 = INTEGER(matDim)[1];
  //is matrix symmetric
  int symmetry = *INTEGER(symmetryR);
  //return sparse form (i.e. list of blocks)?
  int sparse = *INTEGER(sparseR);
  //the index vector that maps locations in the second dimension to
  //first dimension locations, used for non-symmetric case to determine
  //if we need a nugget
  int* ind2_to_1 = INTEGER(ind2_to_1R);
  //extract differentiation flags
  int diff[2];
  extractDiffInfo(diffR, diff);
  
  //loop variables
  int i, j, k, offset, col_offset, n_blocks_tmp;
  //number of parameters (in each block) and flags for where diff happens
  int n_pars;
  int d1 = -1;
  int d2 = -1;
  
  //pointer to list of variables,
  //THIS USES VARIABLE LENGTH ARRAYS, ASSUMING A C99 COMPLIANT COMPILER!!
  double* pars[n_blocks];
  cov_funs_struct cov_funs[n_blocks];
  int* diff_all[n_blocks];
  for(i=0; i<n_blocks; ++i){
    diff_all[i] = Calloc(2, int);
  }

  //return variables
  double *result;
  double *result_all = NULL;
  double* result_part[n_blocks];
  SEXP resultR, result_partR;

  //check length of the inputs
  if( length(parsR)!=n_blocks || length(nuggetR)!=n_blocks ){
    error("Parameters and nugget needs to be of length %d", n_blocks);
  }
  //check if dist is symmetric, when claming to be so.
  if( symmetry!=0 && n1!=n2 ){
    error("Setting the symmetry flag requires a square distance matrix.");
  }
  if( symmetry==0 && length(ind2_to_1R)!=n2){
    //only care about ind2_to_1 in non symmetric case
    error("length(ind2_to_1)!=dim(dist)[2]");
  }
  //extract covariance functions, pointers to parameters and number of parameters
  for(i=0; i<n_blocks; ++i){
    n_pars = length(VECTOR_ELT(parsR, i));
    //check if diff is for this block (diff>0 since we want to exclude 0)
    if( diff[0]>0 && diff[0]<=n_pars ){
      diff_all[i][0] = diff[0];
      d1=i;
    }
    if( diff[1]>0 && diff[1]<=n_pars ){
      diff_all[i][1] = diff[1];
      d2=i;
    }
    //extract covariance function and parameters
    cov_funs[i] = selectCov(CHAR( STRING_ELT(cov_typeR,i) ));
    pars[i] = REAL( VECTOR_ELT(parsR, i) );
    //check that we're getting the right number of parameters
    checkCovarianceFunAndPars(cov_funs[i], CHAR(STRING_ELT(cov_typeR,i)),
                              n_pars, diff_all[i] );
    //adjust diff for possibly matching next block(s)
    diff[0] = diff[0]-n_pars;
    diff[1] = diff[1]-n_pars;
  }//for(i=0; i<n_blocks; ++i)

  if( sparse!=0 ){
    //allocate memory for return list of matrices
    PROTECT(resultR = allocVector(VECSXP, n_blocks));
    for(i=0; i<n_blocks; ++i){
      PROTECT(result_partR = allocMatrix(REALSXP, n1, n2));
      SET_VECTOR_ELT(resultR, i, result_partR);
      result_part[i] = REAL(result_partR);
      //initialise to zero
      memset(result_part[i], 0, n1*n2*sizeof(double) );
    }
    //colum offset is zero since we only do one block at a time
    col_offset = 0;
    n_blocks_tmp = 1;
  }else{
    //allocate memory for return matrix
    PROTECT(resultR = allocMatrix(REALSXP, n1*n_blocks, n2*n_blocks));
    result_all = REAL(resultR);
    //initialise to zero
    memset(result_all, 0, n1*n2*n_blocks*n_blocks*sizeof(double) );
    //colum offset
    col_offset = n1*(n_blocks-1);
    n_blocks_tmp = n_blocks;
  }//if( sparse!=0 ){...}else{...}

  //loop over the blocks, if d1!=d2 we have second derivatives wrt different
  //blocks, i.e. just return zero
  for(i=0; i<n_blocks && (d1==d2 || d2==-1); ++i){
    if(d1!=-1 && d1!=i){
      //compute derivatives (d1 != -1) but all blocks except d1==i are non-zero.
      continue;
    }
    if( sparse!=0 ){
      result = result_part[i];
      offset = 0;
    }else{
      result = result_all;
      offset = i*(n1*n2*n_blocks + n1);
    }
    computeCovarianceBlock(cov_funs[i], pars[i], 0, dist, result + offset,
                           n1, n2, col_offset, symmetry, diff_all[i]);
    //add nugget (if non-zero, and no derivatives)
    if( nugget[i]!=0 && d1==-1){
      if( symmetry!=0 ){
        //symmetric add to all the diagonals
        for(j=0; j<n1; ++j){
          result[offset + j*n1*n_blocks_tmp + j] += nugget[i];
        }
      }else{
        //non-symmetric, use ind2_to_1 to determine nugget
        for(k=0; k<n2; ++k){
          for(j=0; j<n1; ++j){
            if( j==(ind2_to_1[k]-1) ){
              result[offset + k*n1*n_blocks_tmp + j] += nugget[i];
            }
          }
        }
      }//if( symmetry!=0 ){...}else{...}
    }//if( nugget[i]!=0 )
  }//for(i=0; i<n_blocks; ++i)

  //Free the memory allocated for diff
  for(i=0; i<n_blocks; ++i){
    Free( diff_all[i] );
  }
  
  if( sparse!=0 ){
    UNPROTECT(1+n_blocks); /* resultR + list components*/
  }else{
    UNPROTECT(1); /* resultR */
  }
  return resultR;
}//C_makeSigmaB(...)

//compute block cross-covariance matrix for Sigma_nu
SEXP C_makeSigmaNu(SEXP cov_typeR, SEXP parsR, SEXP nuggetR, SEXP randomR,
                   SEXP distR, SEXP blocks1R, SEXP blocks2R, SEXP ind1R,
                   SEXP ind2R, SEXP ind2_to_1R, SEXP symmetryR,
                   SEXP sparseR, SEXP diffR){
  //EXTRACT PARAMETERS FROM R CALL
  //get covariance function name from R variable
  const char* cov_type = CHAR( STRING_ELT(cov_typeR,0) );
  //figure out the function pointers
  cov_funs_struct cov_funs = selectCov(cov_type);
  //pointer to parameters and distance matrix
  double *pars = REAL(parsR);
  double *nugget = REAL(nuggetR);
  int n_nugget = length(nuggetR);
  double random = *REAL(randomR);
  double *dist = REAL(distR);
  int *blocks1 = INTEGER(blocks1R);
  int *blocks2 = INTEGER(blocks2R);
  //dimensions of the distance matrix
  SEXP distDim = getAttrib(distR, R_DimSymbol);
  int n1 = INTEGER(distDim)[0];
  int n2 = INTEGER(distDim)[1];
  //dimensions of the blocks
  int n_blocks = length(blocks1R);
  //the index vectors for locations, used for non-symmetric
  //case to determine if we need a nugget
  int* ind1 = INTEGER(ind1R);
  int* ind2 = INTEGER(ind2R);
  //length of the index vectors, also the size of the final matrix.
  int N1 =length(ind1R);
  int N2 =length(ind2R);
  //ind2_to_1 ignored if symmetry!=0
  int* ind2_to_1 = INTEGER(ind2_to_1R);
  //is matrix symmetric
  int symmetry = *INTEGER(symmetryR);
  //return sparse form (i.e. list of blocks)?
  int sparse = *INTEGER(sparseR);
  //extract differentiation flags
  int diff[2];
  extractDiffInfo(diffR, diff);
  
  //loop variables
  int i, j, i_block;
  //offset in the two dimensions in the current matrix
  int sum1, sum2;
  //sizes of the current block
  int n1_block, n2_block;
  
  //return variables
  double *result = NULL;
  double *covf;
  SEXP resultR, result_partR;

  //check dimensions of cov_type, random_effect and nugget
  if( length(cov_typeR)!=1 || length(randomR)!=1 ){
    error("Sigma_nu assumes same covariance function for all blocks\nlength(cov_type)!=1 || length(random.effect)!=1");
  }
  if( n_nugget!=1 && n_nugget!=n1 ){
    error("Nugget needs to contain 1 or dim(dist)[1] elements.");
  }
  //check if dist is symmetric, when claming to be so.
  if( symmetry!=0 && n1!=n2){
    error("Setting the symmetry flag requires a square distance matrix.");
  }
  //check if covariance function is valid
  checkCovarianceFunAndPars(cov_funs, cov_type, length(parsR), diff);
  //check dimensions of blocks
  if( n_blocks!=length(blocks2R) ){
    error("length(blocks1)!=length(blocks2)");
  }
  sum1 = 0;
  sum2 = 0;
  for(i=0; i<n_blocks; ++i){
    sum1 += blocks1[i];
    sum2 += blocks2[i];
  }
  if( sum1!=N1 ){
    error("sum(blocks1)!=length(ind1)");
  }
  if( sum2!=N2 ){
    error("sum(blocks2)!=length(ind2)");
  }
  //check the index vectors
  if( max_create_cov(ind1, N1)>n1 ){
    error("max(ind1)>dim(dist)[1]");
  }
  if( max_create_cov(ind2, N2)>n2 ){
    error("max(ind2)>dim(dist)[2]");
  }
  if( symmetry==0 && length(ind2_to_1R)!=n2){
    //only care about ind2_to_1 in non symmetric case
    error("length(ind2_to_1)!=dim(dist)[2]");
  }

  //Compute the (cross-)covariance matrix for all locations.
  covf = Calloc(n1*n2, double);
  computeCovarianceBlock(cov_funs, pars, random, dist,
                         covf, n1, n2, 0, symmetry, diff);
  //add nugget (whens length!=1 or nugget!=0, and no diff)
  if( (n_nugget!=1 || nugget[0]!=0) && diff[0]==0){
    if( symmetry!=0 ){
      //symmetric add to the diagonal
      for(i=0; i<n1; ++i){
        if( n_nugget==1 ){
          covf[i*n1 + i] += nugget[0];
        }else{
          covf[i*n1 + i] += nugget[i];
        }
      }//for(i=0; i<n1; ++i){
    }else{
      //non-symmetric, use ind2_to_1 to determine nugget
      for(i=0; i<n2; ++i){
        for(j=0; j<n1; ++j){
          if( j==(ind2_to_1[i]-1) ){
            if( n_nugget==1 ){
              covf[i*n1 + j] += nugget[0];
            }else{
              covf[i*n1 + j] += nugget[j];
            }
          }
        }//for(j=0; j<n1; ++j)
      }//for(i=0; i<n2; ++i)
    }//if( symmetry!=0 ){...}else{...}
  }//if( !(n_nugget!=1 && nugget[0]==0) ){

  //allocate memory for return matrix
  if( sparse!=0 ){
    //allocate memory for return list of matrices
    PROTECT(resultR = allocVector(VECSXP, n_blocks));
  }else{
    PROTECT(resultR = allocMatrix(REALSXP, N1, N2));
    result = REAL(resultR);
    memset(result, 0, N1*N2*sizeof(double));
  }

  //loop over the blocks
  for(i_block=0; i_block<n_blocks; ++i_block){
    //size of the current block
    n1_block = blocks1[i_block];
    n2_block = blocks2[i_block];
    
    if( sparse!=0 ){
      //allocate matrix block for this part
      PROTECT(result_partR = allocMatrix(REALSXP, n1_block, n2_block));
      SET_VECTOR_ELT(resultR, i_block, result_partR);
      result = REAL(result_partR);
      memset(result, 0, n1_block*n2_block*sizeof(double));
      N1 = n1_block;
    }
    
    for(i=0; i<n1_block; ++i){
      for(j=0; j<n2_block; ++j){
        //add relevant element to the large covariance matrix,
        //recall that ind1 and ind2 are base 1 indecies!
        result[i + j*N1] = covf[(ind1[i]-1) + (ind2[j]-1)*n1];
      }
    }
    //incremeant pointers to account for the i:th block
    ind1 += n1_block;
    ind2 += n2_block;
    if( sparse==0 ){
      result += n1_block + N1*n2_block;
    }
  }
  
  //free temporary covariance matrix
  Free(covf);
  if( sparse!=0 ){
    UNPROTECT(1+n_blocks); /* resultR + list components*/
  }else{
    UNPROTECT(1); /* resultR */
  }
  return resultR;
}//C_makeSigmaNu(...)
