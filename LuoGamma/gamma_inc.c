#include "gamma_inc.h"

/* The dominant part,
 * D(a,x) := x^a e^(-x) / Gamma(a+1)
 */
inline static double gamma_inc_D(const double a, const double x, const int index)
{
    double lnr;
    lnr = a * log(x) - x - lngamma_lookup[index]; /*log(gamma(a + 1))*/
    return exp(lnr);
}

/* Continued fraction which occurs in evaluation
 * of Q(a,x) or Gamma(a,x).
 *
 *              1   (1-a)/x  1/x  (2-a)/x   2/x  (3-a)/x
 *   F(a,x) =  ---- ------- ----- -------- ----- -------- ...
 *             1 +   1 +     1 +   1 +      1 +   1 +
 *
 * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no).
 *
 * Split out from gamma_inc_Q_CF() by GJ [Tue Apr  1 13:16:41 MST 2003].
 * See gamma_inc_Q_CF() below.
 *
 */
static double gamma_inc_F_CF(const double a, const double x)
{
    const int    nmax  =  100; /* default 5000 */
    const double small =  DBL_EPSILON*DBL_EPSILON*DBL_EPSILON;

    double hn = 1.0;           /* convergent */
    double Cn = 1.0 / small;
    double Dn = 1.0;
    int n;

    /* n == 1 has a_1, b_1, b_0 independent of a,x,
       so that has been done by hand                */
    for ( n = 2 ; n < nmax ; n++ )
    {
        double an;
        double delta;

        if(GSL_IS_ODD(n))
            an = 0.5*(n-1)/x;
        else
            an = (0.5*n-a)/x;

        Dn = 1.0 + an * Dn;
        if ( fabs(Dn) < small )
            Dn = small;
        Cn = 1.0 + an/Cn;
        if ( fabs(Cn) < small )
            Cn = small;
        Dn = 1.0 / Dn;
        delta = Cn * Dn;
        hn *= delta;

        if(fabs(delta-1.0) < DBL_EPSILON)
        {
            /*printf("%i - continued fraction\n", n);*/

            return hn;
        }
    }

    fprintf(stderr, "\nReached maximum number of iterations a:%e x:%e in %s:%i\n",  a, x, __FILE__, __LINE__);
    abort();

}

/* Continued fraction for Q.
 *
 * Q(a,x) = D(a,x) a/x F(a,x)
 *
 * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no):
 *
 * Since the Gautschi equivalent series method for CF evaluation may lead
 * to singularities, I have replaced it with the modified Lentz algorithm
 * given in
 *
 * I J Thompson and A R Barnett
 * Coulomb and Bessel Functions of Complex Arguments and Order
 * J Computational Physics 64:490-509 (1986)
 *
 * In consequence, gamma_inc_Q_CF_protected() is now obsolete and has been
 * removed.
 *
 * Identification of terms between the above equation for F(a, x) and
 * the first equation in the appendix of Thompson&Barnett is as follows:
 *
 *    b_0 = 0, b_n = 1 for all n > 0
 *
 *    a_1 = 1
 *    a_n = (n/2-a)/x    for n even
 *    a_n = (n-1)/(2x)   for n odd
 *
 */
inline static double gamma_inc_Q_CF(const double a, const double x, const int index)
{
    const double res_D = gamma_inc_D(a, x, index);
    const double res_F = gamma_inc_F_CF(a, x);

    return res_D * (a/x) * res_F;
}


/* P series representation.
 */
static double gamma_inc_P_series(const double a, const double x, const int index)
{
    const int nmax = 100; /* default 10000 */

    double res = gamma_inc_D(a, x, index);

    /* Approximating the terms of the series using Stirling's
       approximation gives t_n = (x/a)^n * exp(-n(n+1)/(2a)), so the
       convergence condition is n^2 / (2a) + (1-(x/a) + (1/2a)) n >>
       -log(GSL_DBL_EPS) if we want t_n < O(1e-16) t_0. The condition
       below detects cases where the minimum value of n is > 5000 */

    /* Normal case: sum the series */


    double sum  = 1.0;
    double term = 1.0;
    int n;

    /* Handle lower part of the series where t_n is increasing, |x| > a+n */

    int nlow = (x > a) ? (x - a): 0;

    for(n=1; n < nlow; n++)
    {
        term *= x/(a+n);
        sum  += term;
    }

    /* Handle upper part of the series where t_n is decreasing, |x| < a+n */

    for (/* n = previous n */ ; n<nmax; n++)
    {
        term *= x/(a+n);
        sum  += term;
        if(fabs(term/sum) < DBL_EPSILON)
        {
            /*printf("%i - Stirling\n", n);*/

            return res * sum;
        }
    }

    fprintf(stderr, "\nReached maximum number of iterations a:%e x:%e in %s:%i\n",  a, x, __FILE__, __LINE__);
    abort();
}

/* valid range x > 0 && x < 1e6 && a > 0 && a < 10 */
double gamma_inc(const double a, const double x, const int index)
{
    const double gamma_a = gamma_lookup[index]; /*gamma(a)*/

    if(x <= 7.) /* magic constant :-) */
    {
        /* If the series is quick, do that. It is
         * robust and simple.
         */

        return gamma_a*(1.0 - gamma_inc_P_series(a, x, index));

    }
    else
    {
        /* Continued fraction is excellent for x >~ a.
               * We do not let x be too large when x > a since
               * it is somewhat pointless to try this there;
               * the function is rapidly decreasing for
               * x large and x > a, and it will just
               * underflow in that region anyway. We
               * catch that case in the standard
               * large-x method.
               */
        return gamma_a*gamma_inc_Q_CF(a, x, index);
    }

}
