// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif


RcppExport SEXP _rcpp_module_boot_stan_fit4BAM_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4FBAM_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4FBAM_MULTI_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4FBAM_MULTI_NF_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4HBAM_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4HBAM_MINI_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4HBAM_MULTI_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4HBAM_MULTI_NF_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4HBAM_NF_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4HBAM_R_MINI_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_stan_fit4BAM_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4BAM_mod, 0},
    {"_rcpp_module_boot_stan_fit4FBAM_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4FBAM_mod, 0},
    {"_rcpp_module_boot_stan_fit4FBAM_MULTI_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4FBAM_MULTI_mod, 0},
    {"_rcpp_module_boot_stan_fit4FBAM_MULTI_NF_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4FBAM_MULTI_NF_mod, 0},
    {"_rcpp_module_boot_stan_fit4HBAM_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4HBAM_mod, 0},
    {"_rcpp_module_boot_stan_fit4HBAM_MINI_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4HBAM_MINI_mod, 0},
    {"_rcpp_module_boot_stan_fit4HBAM_MULTI_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4HBAM_MULTI_mod, 0},
    {"_rcpp_module_boot_stan_fit4HBAM_MULTI_NF_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4HBAM_MULTI_NF_mod, 0},
    {"_rcpp_module_boot_stan_fit4HBAM_NF_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4HBAM_NF_mod, 0},
    {"_rcpp_module_boot_stan_fit4HBAM_R_MINI_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4HBAM_R_MINI_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_hbamr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
