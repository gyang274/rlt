#ifndef __MATRIXMINMAX__
#define __MATRIXMINMAX__


#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// using namespace arma;


//' col_min_idx: colvec min value index within limited range
//' @param wlmt: a limit search on colvec on indices within wlmt
//' @return return an index of min value w.r.t to original index
//' @note cpp use index start from 0 vs r use index start from 1
//' @note in case of equal std:min/std:max take first index seen
// [[Rcpp::export]]
arma::uword col_min_idx(const arma::colvec& u, const arma::ucolvec& wlmt) {

  arma::uword min_val_idx;

  // wlmt.size() == 0 ? u.min(min_val_idx) : ( u(wlmt).min(min_val_idx); min_val_idx = wlmt(min_val_idx); );

  if ( wlmt.size() == 0 ) {

    u.min(min_val_idx);

  } else {

    u(wlmt).min(min_val_idx); min_val_idx = wlmt(min_val_idx);

  }

  return min_val_idx;
}

//' col_max_idx: colvec max value index within limited range
//' @param wlmt: a limit search on colvec on indices within wlmt
//' @return return an index of max value w.r.t to original index
//' @note cpp use index start from 0 vs r use index start from 1
//' @note in case of equal std:min/std:max take first index seen
// [[Rcpp::export]]
arma::uword col_max_idx(const arma::colvec& u, const arma::ucolvec& wlmt) {

  arma::uword max_val_idx;

  // wlmt.size() == 0 ? u.max(max_val_idx);: ( u(wlmt).max(max_val_idx); max_val_idx = wlmt(max_val_idx); );

  if ( wlmt.size() == 0 ) {

    u.max(max_val_idx);

  } else {

    u(wlmt).max(max_val_idx); max_val_idx = wlmt(max_val_idx);

  }

  return max_val_idx;
}

//' col_min_val: colvec min value within limited range
// [[Rcpp::export]]
double col_min_val(const arma::colvec& u, const arma::ucolvec& wlmt) {

  return wlmt.size() > 0 ? u(wlmt).min() : u.min() ;

}

//' col_max_val: colvec max value within limited range
// [[Rcpp::export]]
double col_max_val(const arma::colvec& u, const arma::ucolvec& wlmt) {

  return wlmt.size() > 0 ? u(wlmt).max() : u.max() ;
}

//' col_rgn_val: colvec range = max - min value within limited range
// [[Rcpp::export]]
double col_rgn_val(const arma::colvec& u, const arma::ucolvec& wlmt) {

  return wlmt.size() > 0 ? u(wlmt).max() - u(wlmt).min() : u.max() - u.min() ;
}


#endif // __MATRIXMINMAX__
