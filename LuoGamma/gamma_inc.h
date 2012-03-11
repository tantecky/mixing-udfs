#ifndef GAMMA_INC_H_INCLUDED
#define GAMMA_INC_H_INCLUDED

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

#define GSL_IS_ODD(n)  ((n) & 1)

/*gamma(a)*/
static const double gamma_lookup[] =
    {
        5.07760783719011, /* gamma(2/11)*/
        1.94837447299045, /* gamma(5/11) */
        1.25687274180037 /* gamma(8/11) */
    };

/*log(gamma(a + 1))*/
static const double lngamma_lookup[] =
    {
        -0.0799078397463302, /* log(gamma(2/11 + 1)) */
        -0.121461939023940, /* log(gamma(5/11 + 1)) */
        -0.0898270462533453 /* log(gamma(8/11 + 1)) */
    };



double gamma_inc(const double a, const double x, const int index);

#endif
