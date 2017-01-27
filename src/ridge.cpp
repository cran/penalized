#include <RcppArmadillo.h>
#include <math.h>
#include "core.h"
#include "coxfit.h"
//#include "debug.h"

/* SLOWER THAN R - KEEPING CODE FOR POSSIBLE FUTURE USE */
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List Ridge(   arma::vec beta, arma::vec eta,
                    const arma::mat& Lambda, const arma::mat& X,
                    const Rcpp::Function& fit, const bool trace,
                    const double epsilon, const double maxiter, const Rcpp::List& fitInput )
{
    //CoxFitInput coxIn(fitInput);
    Rcpp::List localfit = Rcpp::as<Rcpp::List>(fit(eta));
    double loglik;
    Rcpp::List W = Rcpp::as<Rcpp::List>(localfit["W"]);
    arma::mat WP = Rcpp::as<arma::mat>(W["P"]);
    arma::vec diagW = Rcpp::as<arma::vec>(W["diagW"]);
    arma::vec residuals = Rcpp::as<arma::vec>(localfit["residuals"]);
    // arma::mat WP;
    // arma::rowvec diagW;
    // arma::rowvec residuals;
    // arma::rowvec fitted;
    // _CoxFitCpp(coxIn, eta.t(), residuals, loglik, WP, diagW, fitted);

    int iter = 0;
    double LL;
    double oldLL = R_NegInf;
    double penalty;
    bool finished = false;

    while (!finished)
    {
        ++iter;

        if (trace)
            Rcpp::Rcout << iter << '\r' << std::flush;

        arma::vec grad;
        if (Lambda.n_cols > 1)
        {
            grad = (X.t() * residuals) - Lambda % beta;
        }
        else
        {
            grad = (X.t() * residuals) - Lambda % beta; // TODO: Test
        }

        arma::mat hessian;
        if (diagW.n_elem == X.n_rows)
        {
            arma::mat XdW(X.n_rows, X.n_cols);
            XdW.each_col() = arma::sqrt(diagW);
            XdW = X % XdW;

            if (WP.n_elem > 1)
            {
                arma::mat wpx = (WP.t() * X);
                hessian =  -(XdW.t() * XdW) + (wpx.t() * wpx);
            }
            else
            {
                hessian = -(XdW.t() * XdW);
            }
        }
        else
        {
            hessian = -(X.t() * X);
        }

        if (Lambda.n_cols > 1)
        {
            hessian = hessian - Lambda;
        }
        else
        {
            hessian.diag() = hessian.diag() - Lambda;
        }

        arma::mat shg = SolveCpp( hessian, grad );
        beta -= shg;
        eta = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X * beta));
        eta = X * beta;

        localfit = Rcpp::as<Rcpp::List>(fit(eta));
        loglik = Rcpp::as<double>(localfit["loglik"]);
        W = Rcpp::as<Rcpp::List>(localfit["W"]);
        WP = Rcpp::as<arma::mat>(W["P"]);
        diagW = Rcpp::as<arma::vec>(W["diagW"]);
        residuals = Rcpp::as<arma::vec>(localfit["residuals"]);
        // _CoxFitCpp(coxIn, eta.t(), residuals, loglik, WP, diagW, fitted);

        if ( Rcpp::NumericVector::is_na(loglik) || iter == maxiter )
        {
            if (trace)
                Rcpp::Rcerr << "Model does not converge: please increase lambda.\n";
            break;
        }

        if (Lambda.n_cols > 1)
        {
            penalty = 0.5 * Sum(beta % (Lambda * beta));
        }
        else
        {
            penalty = 0.5 * Sum(Lambda % arma::square(beta));
        }
        LL = loglik - penalty;

        // Check convergence
        finished = (2 * std::fabs(LL - oldLL) / (2 * std::fabs(LL) + 0.1) < epsilon);
        // half step
        if (LL < oldLL)
        {
            beta = beta + 0.5 * shg;
            eta = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X * beta));
            eta = X * beta;

            localfit = Rcpp::as<Rcpp::List>(fit(eta));
            loglik = Rcpp::as<double>(localfit["loglik"]);
            W = Rcpp::as<Rcpp::List>(localfit["W"]);
            WP = Rcpp::as<arma::mat>(W["P"]);
            diagW = Rcpp::as<arma::vec>(W["diagW"]);
            residuals = Rcpp::as<arma::vec>(localfit["residuals"]);
            // _CoxFitCpp(coxIn, eta.t(), residuals, loglik, WP, diagW, fitted);

            if (Lambda.n_cols > 1)
            {
                penalty = 0.5 * Sum(beta % (Lambda * beta));
            }
            else
            {
                penalty = 0.5 * Sum(Lambda % arma::square(beta));
            }
            LL = loglik - penalty;
        }
        oldLL = LL;
    }

    // Rcpp::List localfit = Rcpp::List::create(   Rcpp::Named("residuals") = residuals,
    //                                             Rcpp::Named("loglik") = loglik,
    //                                             Rcpp::Named("W") = Rcpp::List::create( Rcpp::Named("P") = WP, Rcpp::Named("diagW") = diagW ),
    //                                             Rcpp::Named("lp") = eta.t(),
    //                                             Rcpp::Named("fitted") = fitted
    //                                         );


    return Rcpp::List::create(  Rcpp::Named("beta") = beta,
                                Rcpp::Named("fit") = localfit,
                                Rcpp::Named("penalty") = Rcpp::NumericVector::create( Rcpp::Named("L1") = 0,
                                                                                      Rcpp::Named("L2") = penalty ),
                                Rcpp::Named("iterations") = iter,
                                Rcpp::Named("converged") = finished
                             );
}
