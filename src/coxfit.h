#pragma once
#include <RcppArmadillo.h>

// This struct wraps input needed for performing a cox fit. TODO: Generalize for all models
struct CoxFitInput {
public:
    // Take list given by coxfit.fitInput()
    CoxFitInput( const Rcpp::List& fitInput, const arma::uvec& leftout )
    {
        arma::uvec nLeftoutElem = arma::find(leftout == 0);
        arma::uvec dleftin = arma::find(leftout.elem(Rcpp::as<arma::uvec>(fitInput["whichd"])) == 0);
        status = Rcpp::as<arma::irowvec>(fitInput["status"]).elem(nLeftoutElem);
        riskset = Rcpp::as<arma::umat>(fitInput["Riskset"]).submat(nLeftoutElem, dleftin);
    }
    CoxFitInput( const Rcpp::List& fitInput )
    {
        status = Rcpp::as<arma::irowvec>(fitInput["status"]);
        riskset = Rcpp::as<arma::umat>(fitInput["Riskset"]);
    }

    arma::irowvec status;
    arma::umat riskset;
};

void _CoxFitCpp( const CoxFitInput& input, const arma::rowvec& lp, arma::rowvec& residuals,
    double& loglik, arma::mat& P, arma::rowvec& diagW, arma::rowvec& ws );
