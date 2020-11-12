#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP C_block_mult(SEXP, SEXP, SEXP);
extern SEXP C_calc_F_part_X(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_calc_tFX(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_calc_tFXF(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_cov_names();
extern SEXP C_cov_pars(SEXP);
extern SEXP C_cov_simple(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_dist(SEXP, SEXP, SEXP);
extern SEXP C_dotProd(SEXP, SEXP);
extern SEXP C_inv_chol_block(SEXP, SEXP);
extern SEXP C_make_chol_block(SEXP, SEXP);
extern SEXP C_makeSigmaB(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_makeSigmaNu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_norm2(SEXP);
extern SEXP C_solve_tri_block(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_sum_log(SEXP);
extern SEXP C_sum_log_diag(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"C_block_mult",      (DL_FUNC) &C_block_mult,       3},
    {"C_calc_F_part_X",   (DL_FUNC) &C_calc_F_part_X,    4},
    {"C_calc_tFX",        (DL_FUNC) &C_calc_tFX,         4},
    {"C_calc_tFXF",       (DL_FUNC) &C_calc_tFXF,        5},
    {"C_cov_names",       (DL_FUNC) &C_cov_names,        0},
    {"C_cov_pars",        (DL_FUNC) &C_cov_pars,         1},
    {"C_cov_simple",      (DL_FUNC) &C_cov_simple,       4},
    {"C_dist",            (DL_FUNC) &C_dist,             3},
    {"C_dotProd",         (DL_FUNC) &C_dotProd,          2},
    {"C_inv_chol_block",  (DL_FUNC) &C_inv_chol_block,   2},
    {"C_make_chol_block", (DL_FUNC) &C_make_chol_block,  2},
    {"C_makeSigmaB",      (DL_FUNC) &C_makeSigmaB,       8},
    {"C_makeSigmaNu",     (DL_FUNC) &C_makeSigmaNu,     13},
    {"C_norm2",           (DL_FUNC) &C_norm2,            1},
    {"C_solve_tri_block", (DL_FUNC) &C_solve_tri_block,  4},
    {"C_sum_log",         (DL_FUNC) &C_sum_log,          1},
    {"C_sum_log_diag",    (DL_FUNC) &C_sum_log_diag,     1},
    {NULL, NULL, 0}
};

void R_init_SpatioTemporal(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
