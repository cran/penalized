#pragma once

#include <RcppArmadillo.h>

template <typename T>
double _Sum( T const& obj )
{
    long double result = 0;
    for (auto iter = obj.begin(); iter != obj.end(); ++iter)
    {
        result += *iter;
    }
    return result;
}

// Works the same as R's sum (accu() has less precision)
template <typename T>
inline double Sum( T const& obj )
{
    return _Sum(obj.eval());
}


template <typename T1, typename T2>
inline arma::mat Outer( T1 const& A, T2 const& B )
{
    arma::mat prod(A.n_elem, B.n_elem);
    for (size_t r = 0; r < A.n_elem; ++r)
    {
        for (size_t c = 0; c < B.n_elem; ++c)
        {
            prod.at(r,c) = A.at(r) * B.at(c);
        }
    }
    return prod;
};

inline arma::mat SolveCpp( arma::mat& A, arma::mat& B )
{
    Rcpp::Environment ns = Rcpp::Environment::namespace_env( "penalized" );
    Rcpp::Function solve = ns[".solve"];
    return Rcpp::as<arma::mat>( solve( Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(A)),
        Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(B)) ) );
};


inline arma::mat makeP( const arma::mat& X, const arma::vec& lambda2,
    const arma::vec& lambda1, const arma::vec& signbeta )
{
    int n = X.n_rows;
    int p = X.n_cols;

    arma::uvec free2 = (lambda2 == 0);
    bool free1 = arma::all(lambda1 == 0);
    int m = arma::accu(free2);

    arma::mat P;
    if (free1) {
        P = arma::zeros<arma::mat>(n+m, p);
    } else {
        P = arma::zeros<arma::mat>(n+m+1, p);
    }

    // First columns: free variables in case of no L2-penalization
    arma::uvec free2Elem = arma::find(free2 != 0);
    arma::uvec nFree2Elem = arma::find(free2 == 0);
    for (int i = 0; i < m; ++i)
    {
        P(i, free2Elem[i]) = 1;
    }

    // Next n columns: column span of X
    arma::mat flam(n, p-m);
    flam.each_row() = 1 / lambda2.elem(nFree2Elem).t();
    arma::mat newPs = X.cols(nFree2Elem) % flam;
    arma::vec rows = arma::linspace<arma::vec>(m, n+m-1, n);
    for (size_t ri = 0; ri < rows.n_elem; ++ri)
    {
        for (size_t ci = 0; ci < nFree2Elem.n_elem; ++ci)
        {
            P(rows[ri], nFree2Elem[ci]) = newPs(ri,ci);
        }
    }

    // Additional column due to L1-penalization
    if (!free1) {
        arma::vec newPs2 = (lambda1 % signbeta / lambda2).eval().elem(nFree2Elem);
        for (size_t i = 0; i < nFree2Elem.n_elem; ++i)
        {
            P(n+m, nFree2Elem[i]) = newPs2[i];
        }
    }

    arma::uword rank = arma::rank(P);
    if (rank < P.n_rows)
    {
        arma::mat Q, R;
        arma::qr(Q, R, P.t());
        arma::umat keep = arma::sort_index( arma::abs(R.diag()), "descend" ).eval().elem( arma::linspace<arma::uvec>(0, rank-1, rank) );
        P = P.rows( keep );
    }

    return P;
};
