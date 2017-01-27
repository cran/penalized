#include <RcppArmadillo.h>
#include <math.h>
#include "core.h"
#include "coxfit.h"

/*
    Avoids extra copy of beta and does not return anything. Only to be called from within cpp.
*/
void _Lasso( arma::vec& beta, const arma::vec& lambda,
                    const arma::vec& lambda2, const arma::uvec& positive,
                    const arma::mat& X, const Rcpp::Function& fit, const bool trace,
                    const double epsilon, const double maxiter )
{
    arma::vec grad;

    bool enet = arma::any(lambda2 != 0);
    bool finished = false;
    bool finishedLL = false;
    bool finishedpen = false;
    bool newfit = true;
    bool whereNRSet = false;
    bool tryNR = false;
    bool NRfailed = false;

    int m = beta.n_elem;
    int n = X.n_rows;

    double LL = R_NegInf;
    double penalty = R_PosInf;
    double penalty1 = R_PosInf;
    double penalty2 = R_PosInf;
    double retain = 0.05;
    double oldLL = 0;
    double oldPenalty = 0;
    double topt = 0;
    double cumsteps = 0;

    int iter = 0;
    int nvar = m;

    arma::mat activeX;
    arma::mat P;
    arma::mat PX;
    arma::mat PlP;
    arma::mat gams;
    arma::mat NRbeta;

    arma::uvec active = arma::ones<arma::uvec>(m);
    arma::uvec activeElem = arma::find( active );
    arma::uvec _free = ((lambda == 0) + (positive == 0)) == 2;
    arma::uvec whereNR;

    // localfit variables
    Rcpp::List localfit;
    double loglik;
    arma::mat WP;
    arma::rowvec diagW;
    arma::vec residuals;

    if (trace)
        Rcpp::Rcout << "# nonzero coefficients: " << m << '\r' << std::flush;

    while (!finished)
    {
        arma::uvec nzb = beta != 0;
        arma::uvec nzbElem = arma::find(nzb);

        // Calculate local-likelihood fit
        if (newfit)
        {
            activeX = X.cols( nzbElem );

            arma::vec linpred;
            if ( activeX.is_empty() )
            {
                linpred = arma::zeros<arma::vec>(X.n_rows);
            }
            else
            {
                linpred = (activeX * beta.elem( nzbElem ));
            }

            localfit = Rcpp::as<Rcpp::List>(fit( Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(linpred)) ));
            loglik = Rcpp::as<double>(localfit["loglik"]);
            Rcpp::List W = Rcpp::as<Rcpp::List>(localfit["W"]);
            WP = Rcpp::as<arma::mat>(W["P"]);
            diagW = Rcpp::as<arma::rowvec>(W["diagW"]);
            residuals = Rcpp::as<arma::vec>(localfit["residuals"]);
            //_CoxFitCpp(fitInput, linpred, residuals, loglik, WP, diagW);

            // Check for divergence
            if ( Rcpp::NumericVector::is_na(loglik) )
            {
                if (trace)
                    Rcpp::Rcerr << "Model does not converge: please increase lambda.\n"; // TODO: Use R's warning() function
                break;
            }

            grad = (X.t() * residuals);

            if (enet)
            {
                grad.elem( activeElem ) = grad.elem( activeElem ) -
                    lambda2.elem( activeElem ) % beta.elem( activeElem );
            }

            oldLL = LL;
            oldPenalty = penalty;
            LL = loglik;
            penalty1 = Sum( lambda.elem(activeElem) % arma::abs(beta.elem(activeElem)));

            if (enet)
            {
                penalty2 = Sum( lambda2.elem( activeElem ) % arma::square(beta.elem( activeElem )) );
            }
            else
            {
                penalty2 = 0;
            }


            penalty = penalty1 + penalty2;
            finishedLL = (2 * std::fabs(LL - oldLL) / (2 * std::fabs(LL - penalty) + 0.1) < epsilon);
            finishedpen = (2 * std::fabs(penalty - oldPenalty) / (2 * std::fabs(LL - penalty) + 0.1) < epsilon);
            cumsteps = 0;
        }

        // Calculate the penalized gradient from the likelihood gradient
        arma::vec direction(m, arma::fill::zeros);
        direction.elem( nzbElem ) = grad.elem(nzbElem) -
            ( lambda.elem(nzbElem) % arma::sign(beta.elem(nzbElem)) );

        arma::uvec newb(nzb.n_elem);
        for (size_t i = 0; i != newb.n_elem; ++i)
        {
            newb[i] = (arma::uword)( (nzb[i] == 0) && (positive[i] ?
                (grad[i] > lambda[i]) : (std::fabs(grad[i]) > lambda[i])) );
        }

        arma::uvec newbElem = arma::find(newb);
        direction.elem( newbElem ) = grad.elem(newbElem) - lambda.elem(newbElem) % arma::sign(grad.elem(newbElem));
        arma::uvec oldActive = active;
        active = (nzb + newb) != 0;
        activeElem = arma::find( active );
        arma::vec activeBeta = beta.elem( activeElem );
        arma::vec activeDir = direction.elem( activeElem );

        // check if retaining the old fit of the model does more harm than good
        int oldnvar = nvar;
        nvar = Sum(active);
        if ((oldLL - oldPenalty > LL - penalty) || (nvar > 1.1 * oldnvar))
        {
            retain *= 0.5;
        }

        // Check convergence
        bool finishednvar = !arma::any(active != oldActive);
        finished = (finishedLL && finishedpen && finishednvar) || arma::all(activeDir == 0) || (iter == maxiter);

        if (!finished) {
            ++iter;

            if (tryNR) {
                activeX = X.cols( activeElem );

                if ( enet && (nvar > (n + 1 + Sum(_free))) )
                {
                    if ( !whereNRSet || arma::any(whereNR != active) )
                    {
                        whereNRSet = true;
                        whereNR = active;
                        P = makeP(activeX, lambda2.elem(activeElem),
                            lambda.elem(activeElem), arma::sign(activeBeta));
                        arma::mat gamsA = (P * P.t());
                        arma::mat gamsB = (P * activeBeta);
                        gams = SolveCpp( gamsA, gamsB );
                        PX = P * activeX.t();

                        arma::mat Pl(P.n_rows, P.n_cols);
                        Pl.each_row() = arma::sqrt(lambda2.elem(activeElem)).t();
                        Pl = P % Pl;
                        PlP = Pl * Pl.t();
                    }

                    arma::mat hessian;

                    if (diagW.n_elem == PX.n_cols)
                    {
                        arma::rowvec diagWRt = arma::sqrt(diagW);
                        arma::mat PXdW(PX.n_rows, PX.n_cols);
                        PXdW.each_row() = diagWRt;
                        PXdW = PX % PXdW;

                        if (WP.n_elem > 1)
                        {
                            arma::mat ppx = (PX * WP).t();
                            hessian =  -(PXdW * PXdW.t()) + (ppx.t() * ppx) - PlP;
                        }
                        else
                        {
                            hessian = -(PXdW * PXdW.t()) - PlP;
                        }
                    }
                    else
                    {
                        hessian = -(PX * PX.t()) - PlP;
                    }

                    arma::mat Pgrad = P * direction.elem(activeElem);
                    arma::mat shg = SolveCpp( hessian, Pgrad );
                    gams = gams - shg;
                    NRbeta = P.t() * gams;

                }
                else
                {
                    arma::mat hessian;

                    if (diagW.n_elem == activeX.n_rows)
                    {
                        arma::vec diagWRt = arma::sqrt(diagW).t();
                        arma::mat XdW(activeX.n_rows, activeX.n_cols);
                        XdW.each_col() = diagWRt;
                        XdW = (activeX % XdW);

                        if (WP.n_elem > 1)
                        {
                            arma::mat cProdPX = WP.t() * activeX;
                            hessian = -(XdW.t() * XdW) + (cProdPX.t() * cProdPX);
                        }
                        else
                        {
                            hessian = -(XdW.t() * XdW);
                        }
                    }
                    else
                    {
                        hessian = -(activeX.t() * activeX);
                    }

                    if (enet)
                    {
                        hessian.diag() = hessian.diag() - lambda2.elem(activeElem);
                    }

                    NRbeta = activeBeta - SolveCpp(hessian, activeDir);
                }

                NRfailed = !arma::all(arma::sign(NRbeta) == arma::sign(activeBeta));

                if (!NRfailed)
                {
                    beta.elem(activeElem) = NRbeta;
                    newfit = true;
                }
            }

            if (!tryNR || NRfailed)
            {
                // find the second derivative of the likelihood in the projected direction
                if (newfit)
                {
                    arma::mat Xdir = X.cols(activeElem) * activeDir;
                    double curve;

                    if (diagW.n_elem > 0)
                    {
                        if (WP.n_elem > 1)
                        {
                            arma::mat cpPX = WP.t() * Xdir;
                            curve =  (Sum( arma::square(Xdir).t() % diagW ) -
                                arma::as_scalar(cpPX.t() * cpPX)) / Sum( arma::square(activeDir) );
                        }
                        else
                        {
                            curve = Sum( arma::square(Xdir) % diagW ) / Sum( arma::square(activeDir) );
                        }
                    }
                    else
                    {
                        curve = Sum( arma::square(Xdir) ) / Sum( arma::square(activeDir) );
                    }

                    if (enet)
                    {
                        curve = curve + Sum(lambda2.elem(activeElem) %
                            arma::square(activeDir)) / Sum( arma::square(activeDir) );
                    }

                    topt = 1 / curve;
                }

                // how far can we go in the calculated direction before finding a new zero?
                arma::vec tedge(m, arma::fill::zeros);
                tedge.elem( activeElem ) = -activeBeta / activeDir;
                tedge.elem( arma::find(tedge <= 0) ).fill(2 * topt);
                tedge.elem( arma::find(_free > 0) ).fill(2 * topt);
                double mintedge = arma::min(tedge);

                // recalculate beta
                if (mintedge + cumsteps < topt)
                {
                    beta.elem(activeElem) = activeBeta + (mintedge * activeDir);
                    beta.elem( arma::find(tedge == mintedge) ).zeros();
                    cumsteps = cumsteps + mintedge;
                    newfit = (cumsteps > (retain * topt)) || (nvar == 1);
                    NRfailed = false;
                    tryNR = false;
                }
                else
                {
                    beta.elem( activeElem ) = (activeBeta + (topt - cumsteps) * activeDir);
                    tryNR = (cumsteps == 0) && !NRfailed && finishednvar && (enet || nvar < n);
                    newfit = true;
                }
            }
        }

        if (trace)
        {
            Rcpp::Rcout << "# nonzero coefficients: " << nvar << '\r' << std::flush;
        }
    }
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List Lasso( arma::vec beta, const arma::vec& lambda,
                    const arma::vec& lambda2, const arma::uvec& positive,
                    const arma::mat& X, const Rcpp::Function& fit, const bool trace,
                    const double epsilon, const double maxiter )
{
    arma::vec grad;

    bool enet = arma::any(lambda2 != 0);
    bool finished = false;
    bool finishedLL = false;
    bool finishedpen = false;
    bool newfit = true;
    bool converged;
    bool whereNRSet = false;
    bool tryNR = false;
    bool NRfailed = false;

    int m = beta.n_elem;
    int n = X.n_rows;

    double LL = R_NegInf;
    double penalty = R_PosInf;
    double penalty1 = R_PosInf;
    double penalty2 = R_PosInf;
    double retain = 0.05;
    double oldLL = 0;
    double oldPenalty = 0;
    double topt = 0;
    double cumsteps = 0;

    int iter = 0;
    int nvar = m;

    arma::mat activeX;
    arma::mat P;
    arma::mat PX;
    arma::mat PlP;
    arma::mat gams;
    arma::mat NRbeta;

    arma::uvec active = arma::ones<arma::uvec>(m);
    arma::uvec activeElem = arma::find( active );
    arma::uvec _free = ((lambda == 0) + (positive == 0)) == 2;
    arma::uvec whereNR;

    // localfit variables
    Rcpp::List localfit;
    double loglik;
    arma::mat WP;
    arma::rowvec diagW;
    arma::vec residuals;

    if (trace)
        Rcpp::Rcout << "# nonzero coefficients: " << m << '\r' << std::flush;

    while (!finished)
    {
        arma::uvec nzb = beta != 0;
        arma::uvec nzbElem = arma::find(nzb);

        // Calculate local-likelihood fit
        if (newfit)
        {
            activeX = X.cols( nzbElem );

            arma::vec linpred;
            if ( activeX.is_empty() )
            {
                linpred = arma::zeros<arma::vec>(X.n_rows);
            }
            else
            {
                linpred = (activeX * beta.elem( nzbElem ));
            }

            localfit = Rcpp::as<Rcpp::List>(fit( Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(linpred)) ));
            loglik = Rcpp::as<double>(localfit["loglik"]);
            Rcpp::List W = Rcpp::as<Rcpp::List>(localfit["W"]);
            WP = Rcpp::as<arma::mat>(W["P"]);
            diagW = Rcpp::as<arma::rowvec>(W["diagW"]);
            residuals = Rcpp::as<arma::vec>(localfit["residuals"]);

            // Check for divergence
            if ( Rcpp::NumericVector::is_na(loglik) )
            {
                if (trace)
                    Rcpp::Rcerr << "Model does not converge: please increase lambda.\n"; // TODO: Use R's warning() function

                converged = false;
                break;
            }

            grad = (X.t() * residuals);

            if (enet)
            {
                grad.elem( activeElem ) = grad.elem( activeElem ) -
                    lambda2.elem( activeElem ) % beta.elem( activeElem );
            }

            oldLL = LL;
            oldPenalty = penalty;
            LL = loglik;
            penalty1 = Sum( lambda.elem(activeElem) % arma::abs(beta.elem(activeElem)));

            if (enet)
            {
                penalty2 = Sum( lambda2.elem( activeElem ) % arma::square(beta.elem( activeElem )) );
            }
            else
            {
                penalty2 = 0;
            }

            penalty = penalty1 + penalty2;

            finishedLL = (2 * std::fabs(LL - oldLL) / (2 * std::fabs(LL - penalty) + 0.1) < epsilon);
            finishedpen = (2 * std::fabs(penalty - oldPenalty) / (2 * std::fabs(LL - penalty) + 0.1) < epsilon);
            cumsteps = 0;
        }

        // Calculate the penalized gradient from the likelihood gradient
        arma::vec direction(m, arma::fill::zeros);
        direction.elem( nzbElem ) = grad.elem(nzbElem) -
            ( lambda.elem(nzbElem) % arma::sign(beta.elem(nzbElem)) );

        arma::uvec newb(nzb.n_elem);
        for (size_t i = 0; i != newb.n_elem; ++i)
        {
            newb[i] = (arma::uword)( (nzb[i] == 0) && (positive[i] ?
                (grad[i] > lambda[i]) : (std::fabs(grad[i]) > lambda[i])) );
        }

        arma::uvec newbElem = arma::find(newb);
        direction.elem( newbElem ) = grad.elem(newbElem) - lambda.elem(newbElem) % arma::sign(grad.elem(newbElem));
        arma::uvec oldActive = active;
        active = (nzb + newb) != 0;
        activeElem = arma::find( active );
        arma::vec activeBeta = beta.elem( activeElem );
        arma::vec activeDir = direction.elem( activeElem );

        // check if retaining the old fit of the model does more harm than good
        int oldnvar = nvar;
        nvar = Sum(active);
        if ((oldLL - oldPenalty > LL - penalty) || (nvar > 1.1 * oldnvar))
        {
            retain *= 0.5;
        }

        // Check convergence
        bool finishednvar = !arma::any(active != oldActive);
        finished = (finishedLL && finishedpen && finishednvar) || arma::all(activeDir == 0) || (iter == maxiter);

        if (!finished) {
            ++iter;

            if (tryNR) {
                activeX = X.cols( activeElem );

                if ( enet && (nvar > (n + 1 + Sum(_free))) )
                {
                    if ( !whereNRSet || arma::any(whereNR != active) )
                    {
                        whereNRSet = true;
                        whereNR = active;
                        P = makeP(activeX, lambda2.elem(activeElem),
                            lambda.elem(activeElem), arma::sign(activeBeta));
                        arma::mat gamsA = (P * P.t());
                        arma::mat gamsB = (P * activeBeta);
                        gams = SolveCpp( gamsA, gamsB );
                        PX = P * activeX.t();

                        arma::mat Pl(P.n_rows, P.n_cols);
                        Pl.each_row() = arma::sqrt(lambda2.elem(activeElem)).t();
                        Pl = P % Pl;
                        PlP = Pl * Pl.t();
                    }

                    arma::mat hessian;

                    if (diagW.n_elem == PX.n_cols)
                    {
                        arma::rowvec diagWRt = arma::sqrt(diagW);
                        arma::mat PXdW(PX.n_rows, PX.n_cols);
                        PXdW.each_row() = diagWRt;
                        PXdW = PX % PXdW;

                        if (WP.n_elem > 1)
                        {
                            arma::mat ppx = (PX * WP).t();
                            hessian =  -(PXdW * PXdW.t()) + (ppx.t() * ppx) - PlP;
                        }
                        else
                        {
                            hessian = -(PXdW * PXdW.t()) - PlP;
                        }
                    }
                    else
                    {
                        hessian = -(PX * PX.t()) - PlP;
                    }
                    arma::mat Pgrad = P * direction.elem(activeElem);
                    arma::mat shg = SolveCpp( hessian, Pgrad );
                    gams = gams - shg;
                    NRbeta = P.t() * gams;

                }
                else
                {
                    arma::mat hessian;

                    if (diagW.n_elem == activeX.n_rows)
                    {
                        arma::vec diagWRt = arma::sqrt(diagW).t();
                        arma::mat XdW(activeX.n_rows, activeX.n_cols);
                        XdW.each_col() = diagWRt;
                        XdW = (activeX % XdW);

                        if (WP.n_elem > 1)
                        {
                            arma::mat cProdPX = WP.t() * activeX;
                            hessian = -(XdW.t() * XdW) + (cProdPX.t() * cProdPX);
                        }
                        else
                        {
                            hessian = -(XdW.t() * XdW);
                        }
                    }
                    else
                    {
                        hessian = -(activeX.t() * activeX);
                    }

                    if (enet)
                    {
                        hessian.diag() = hessian.diag() - lambda2.elem(activeElem);
                    }

                    NRbeta = activeBeta - SolveCpp(hessian, activeDir);
                }

                NRfailed = !arma::all(arma::sign(NRbeta) == arma::sign(activeBeta));

                if (!NRfailed)
                {
                    beta.elem(activeElem) = NRbeta;
                    newfit = true;
                }
            }

            if (!tryNR || NRfailed)
            {
                // find the second derivative of the likelihood in the projected direction
                if (newfit)
                {
                    arma::mat Xdir = X.cols(activeElem) * activeDir;
                    double curve;

                    if (diagW.n_elem > 0)
                    {
                        if (WP.n_elem > 1)
                        {
                            arma::mat cpPX = WP.t() * Xdir;

                            curve =  (Sum( arma::square(Xdir).t() % diagW ) -
                                (cpPX.t() * cpPX).eval()[0]) / Sum( arma::square(activeDir) ); // as_scalar replaced for accuracy
                        }
                        else
                        {
                            curve = Sum( arma::square(Xdir) % diagW.t() ) / Sum( arma::square(activeDir) );
                        }
                    }
                    else
                    {
                        curve = Sum( arma::square(Xdir) ) / Sum( arma::square(activeDir) );
                    }

                    if (enet)
                    {
                        curve = curve + Sum(lambda2.elem(activeElem) %
                            arma::square(activeDir)) / Sum( arma::square(activeDir) );
                    }

                    topt = 1 / curve;
                }

                // how far can we go in the calculated direction before finding a new zero?
                arma::vec tedge(m, arma::fill::zeros);
                tedge.elem( activeElem ) = -activeBeta / activeDir;
                tedge.elem( arma::find(tedge <= 0) ).fill(2 * topt);
                tedge.elem( arma::find(_free > 0) ).fill(2 * topt);
                double mintedge = arma::min(tedge);

                // recalculate beta
                if (mintedge + cumsteps < topt)
                {
                    beta.elem(activeElem) = activeBeta + (mintedge * activeDir);
                    beta.elem( arma::find(tedge == mintedge) ).zeros();
                    cumsteps = cumsteps + mintedge;
                    newfit = (cumsteps > (retain * topt)) || (nvar == 1);
                    NRfailed = false;
                    tryNR = false;
                }
                else
                {
                    beta.elem( activeElem ) = (activeBeta + (topt - cumsteps) * activeDir);
                    tryNR = (cumsteps == 0) && !NRfailed && finishednvar && (enet || nvar < n);
                    newfit = true;
                }
            }
        }
        else
        {
            converged = (iter < maxiter);
        }

        if (trace)
        {

            Rcpp::Rcout << "\r# nonzero coefficients: " << nvar <<  "          " << '\r' << std::flush;
        }
    }
    return Rcpp::List::create(  Rcpp::Named("beta") = beta,
                                Rcpp::Named("fit") = localfit,
                                Rcpp::Named("penalty") = Rcpp::NumericVector::create( Rcpp::Named("L1") = penalty1,
                                                                                      Rcpp::Named("L2") = penalty2 ),
                                Rcpp::Named("iterations") = iter,
                                Rcpp::Named("converged") = converged
                            );
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List StepLasso(   arma::vec beta, const arma::vec& lambda,
                        const arma::vec& lambda2, const arma::uvec& positive,
                        const arma::mat& X, const Rcpp::Function& fit, const bool trace,
                        const double epsilon, const double maxiter )
{
    double nextLambda;
    size_t n = X.n_rows;
    bool finished = false;

    while (true)
    {
        arma::uvec nzbElem = arma::find(beta != 0);
        arma::uvec zbElem = arma::find(beta == 0);

        // Calculate local-likelihood fit
        arma::mat activeX = X.cols( nzbElem );

        arma::vec linpred;
        if ( activeX.is_empty() )
        {
            linpred = arma::zeros<arma::vec>(n);
        }
        else
        {
            linpred = (activeX * beta.elem( nzbElem ));
        }

        Rcpp::List localfit = Rcpp::as<Rcpp::List>(fit( Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(linpred)) ));
        arma::vec residuals = Rcpp::as<arma::vec>(localfit["residuals"]);

        arma::vec grad = (X.cols(zbElem).t() * residuals);
        arma::vec rel = grad / lambda.elem(zbElem);
        rel = rel.elem( arma::find((rel > 0) + (positive.elem(zbElem) == 0) > 0) );

        if (rel.n_elem > n)
        {
            nextLambda = arma::sort(arma::abs(rel), "descend").eval()[n-1];
        }
        else
        {
            finished = true;
        }

        if (!finished && nextLambda > 1)
        {
            _Lasso(beta, (nextLambda * lambda), lambda2, positive, X, fit, trace, 1e-4, maxiter);
        }
        else
        {
            return(Lasso(beta, lambda, lambda2, positive, X, fit, trace, epsilon, maxiter));
        }
    }
}
