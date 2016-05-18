// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ChangePoints_fromR
IntegerVector ChangePoints_fromR(const NumericVector& x);
RcppExport SEXP pvar_ChangePoints_fromR(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    __result = Rcpp::wrap(ChangePoints_fromR(x));
    return __result;
END_RCPP
}
// pvarC
List pvarC(const NumericVector& x, double& p, int LSI);
RcppExport SEXP pvar_pvarC(SEXP xSEXP, SEXP pSEXP, SEXP LSISEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double& >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type LSI(LSISEXP);
    __result = Rcpp::wrap(pvarC(x, p, LSI));
    return __result;
END_RCPP
}
// AddPvar
List AddPvar(List PV1, List PV2, bool AddIfPossible);
RcppExport SEXP pvar_AddPvar(SEXP PV1SEXP, SEXP PV2SEXP, SEXP AddIfPossibleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type PV1(PV1SEXP);
    Rcpp::traits::input_parameter< List >::type PV2(PV2SEXP);
    Rcpp::traits::input_parameter< bool >::type AddIfPossible(AddIfPossibleSEXP);
    __result = Rcpp::wrap(AddPvar(PV1, PV2, AddIfPossible));
    return __result;
END_RCPP
}
// test_CheckSmallIntervals
NumericVector test_CheckSmallIntervals(const NumericVector& x, const double& p, const int& dn);
RcppExport SEXP pvar_test_CheckSmallIntervals(SEXP xSEXP, SEXP pSEXP, SEXP dnSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int& >::type dn(dnSEXP);
    __result = Rcpp::wrap(test_CheckSmallIntervals(x, p, dn));
    return __result;
END_RCPP
}
// test_prepare_prt
List test_prepare_prt(const NumericVector& x, const double& p);
RcppExport SEXP pvar_test_prepare_prt(SEXP xSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type p(pSEXP);
    __result = Rcpp::wrap(test_prepare_prt(x, p));
    return __result;
END_RCPP
}
