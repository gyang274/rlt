#include "matrix_minmax.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// skm_gdp_val0_cpp: selective skm with greedy propgation on val0 cells.
//
// skm_gdp_val0_cpp with an input m x n matrix: matrix is init as all rows indentical
// and sum to 1 followed by mask a selective columns for each row as 0 - m are stores
// sources, n are zipcode desinations and v the cell values are percentage population
// lives in the destination, and cells are masked as 0 when s - t within a given dist
// X miles, want to select m sequentially so that always covers as much as popluation
// objective is sequential max(sum( min(dist(s, t) - all s) < Xmile ? p : 0 - all t))
//
// @param x: m x n matrix: s - t - p
//
// @return s: index 1 - m gives highest overall p:
//   max(sum( min(dist(s, t) - all s) < Xmile ? p : 0 - all t))

// [[Rcpp::export]]
arma::uvec skm_gdp_val0_cpp(const arma::mat& x) {

  // Rcpp::Rcout << "skm_gdp_val0_cpp: init x as:" << std::endl <<  x   << std::endl;

  // init output vector s
  arma::uvec s = arma::zeros<arma::uvec>(x.n_rows);

  // Rcpp::Rcout << "skm_gdp_val0_cpp: init s as:" << std::endl << s.t() << std::endl;

  arma::uvec ulmt = arma::linspace<arma::uvec>(0, x.n_rows - 1, x.n_rows);

  arma::uvec vlmt = arma::linspace<arma::uvec>(0, x.n_cols - 1, x.n_cols);

  // Rcpp::Rcout << "skm_gdp_val0_cpp: init ulmt as:" << std::endl << ulmt.t() << std::endl;

  // Rcpp::Rcout << "skm_gdp_val0_cpp: init vlmt as:" << std::endl << vlmt.t() << std::endl;

  for ( arma::uword i = 0; i < x.n_rows; i++ ) {

    Rcpp::Rcout << "skm_gdp_val0_cpp: optimize at it <" << i << "> ..." << std::endl;

    if ( vlmt.size() > 0 ) {

      // rowSum %p outside Xmile for each location take into account s selected
      arma::vec rs = arma::sum(x.cols(vlmt), 1);

      // Rcpp::Rcout << "skm_gdp_val0_cpp: optimize at it <" << i << "> ... update rs: " << std::endl << rs.t() << std::endl;

      s(i) = col_min_idx(rs, ulmt);

      // Rcpp::Rcout << "skm_gdp_val0_cpp: optimize at it <" << i << "> ... update s: " << std::endl << s.t() << std::endl;

      // update ulmt and vlmt
      ulmt = ulmt(arma::find(ulmt != s(i)));

      // Rcpp::Rcout << "skm_gdp_val0_cpp: optimize at it <" << i << "> ... update ulmt: " << std::endl << ulmt.t() << std::endl;

      arma::vec u = x.row(s(i)).t();

      vlmt = vlmt(arma::find( u(vlmt) != 0 ));

      // Rcpp::Rcout << "skm_gdp_val0_cpp: optimize at it <" << i << "> ... update vlmt: " << std::endl << vlmt.t() << std::endl;

    } else {

      // Rcpp::Rcout << "skm_gdp_val0_cpp: optimize at it <" << i << "> ... break ..." << std::endl;

      break;

    }

  }

  if ( ulmt.size() > 0 ) {

    s.tail(ulmt.size()) = ulmt;

  }

  return s;
}


// skm_gdp_val0_mt_cpp: create metric view of skm_gdp_val0 output s on %pop matrix x
//
// @param x: matrix of cost - row s sum should be 1 and then set t within Xmile to 0
//           column t cost should be the same expect for s within Xmile was set to 0
// @param s: list of s index being selected sequentially into output by skm_gdp_val0
//
// @return d: incremential cost deduction (cost => 0) by adding i-th s into the list

// [[Rcpp::export]]
arma::vec skm_gdp_val0_mt_cpp(const arma::mat& x, const arma::uvec& s) {

  // Rcpp::Rcout << "skm_gdp_val0_mt_cpp: init x as:" << std::endl <<  x   << std::endl;

  // Rcpp::Rcout << "skm_gdp_val0_mt_cpp: init s as:" << std::endl << s.t() << std::endl;

  // init output vector d
  arma::vec d = arma::zeros<arma::vec>(s.size());

  // Rcpp::Rcout << "skm_gdp_val0_mt_cpp: init d as:" << std::endl << d.t() << std::endl;

  arma::uvec vlmt = arma::linspace<arma::uvec>(0, x.n_cols - 1, x.n_cols);

  for (arma::uword i = 0; i < s.size(); i++ ) {

    Rcpp::Rcout << "skm_gdp_val0_mt_cpp: optimize at it <" << i << "> ..." << std::endl;

    // calculate cumulative cost deduction d and update vlmt
    arma::vec u = x.row(s(i)).t();

    d(i) = 1.0 - arma::sum(u(vlmt));

    // Rcpp::Rcout << "skm_gdp_val0_mt_cpp: calculate at i <" << i << "> - update d: " << std::endl << d.t() << std::endl;

    vlmt = vlmt(arma::find( u(vlmt) != 0 ));

    // Rcpp::Rcout << "skm_gdp_val0_cpp: optimize at it <" << i << "> ... update vlmt: " << std::endl << vlmt.t() << std::endl;

  }

  return d;

}

