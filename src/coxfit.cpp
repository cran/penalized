#include <RcppArmadillo.h>
#include "core.h"
#include "coxfit.h"

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List CoxFitCpp(   const arma::rowvec& lp, const arma::irowvec&  status,
                        const arma::umat& riskset )
{
    arma::rowvec ws = arma::exp(lp);

    bool exploded = arma::any(ws == 0 || ws == R_PosInf);

    if (exploded)
    {
        ws.fill(1e+10);
        ws.elem( arma::find(status == 0) ).fill(1e-10);
    }

    arma::urowvec onesVec( riskset.n_cols, arma::fill::ones );
    arma::rowvec breslows = onesVec / (ws * riskset);
    arma::colvec breslow = riskset * breslows.t();

    // The martingale residuals
    arma::rowvec wbres = breslow.t() % ws; // Operator chaining (a_status - breslow.t() % a_ws) causes segfault
    arma::rowvec residuals = status - wbres;

    //logArma(wbres.t());

    // Log-likelihood
    double loglik = NA_REAL;
    if (!exploded)
        loglik = -Sum(wbres) + Sum(arma::log(breslows)) +
                Sum( lp.elem(arma::find(status == 1)) );
    if (loglik == R_NegInf) loglik = NA_REAL; // Needed? What cases?

    // Weights
    arma::mat Pij = Outer(ws, breslows);
    Pij.elem( arma::find(riskset == 0) ).zeros();
    Rcpp::List W = Rcpp::List::create( Rcpp::Named("P") = Pij, Rcpp::Named("diagW") = wbres );

    // Re-adjust ws for output if expolded
    if (exploded)
        ws = arma::exp(lp);

    return Rcpp::List::create(  Rcpp::Named("residuals") = residuals,
                                Rcpp::Named("loglik") = loglik,
                                Rcpp::Named("W") = W,
                                Rcpp::Named("lp") = lp,
                                Rcpp::Named("fitted") = ws
                                );
}


// Full c++ cox fit. All params besides input for output only
void _CoxFitCpp( const CoxFitInput& input, const arma::rowvec& lp, arma::rowvec& residuals,
    double& loglik, arma::mat& P, arma::rowvec& diagW, arma::rowvec& ws )
{
    ws = arma::exp(lp);

    bool exploded = arma::any(ws == 0 || ws == R_PosInf);

    if (exploded)
    {
        ws.fill(1e+10);
        ws.elem( arma::find(input.status == 0) ).fill(1e-10);
    }

    arma::urowvec onesVec( input.riskset.n_cols, arma::fill::ones );
    arma::rowvec breslows = onesVec / (ws * input.riskset);
    arma::colvec breslow = input.riskset * breslows.t();

    // The martingale residuals
    diagW = breslow.t() % ws; // Operator chaining (a_status - breslow.t() % a_ws) causes segfault
    residuals = input.status - diagW;

    // Log-likelihood
    loglik = NA_REAL;
    if (!exploded)
        loglik = -Sum(diagW) + Sum(arma::log(breslows)) +
                Sum( lp.elem(arma::find(input.status == 1)) );
    if (loglik == R_NegInf) loglik = NA_REAL; // Needed? What cases?

    // Weights
    P = Outer(ws, breslows);
    P.elem( arma::find(input.riskset == 0) ).zeros();

    // Re-adjust ws for output if expolded
    if (exploded)
        ws = arma::exp(lp);
}
