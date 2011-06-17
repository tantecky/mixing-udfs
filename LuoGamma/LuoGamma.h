#ifndef LUO_H_INCLUDED
#define LUO_H_INCLUDED

#include "udf.h"
#include "metric.h"
#include "sg_pb.h"
#include "sg_mphase.h"
#include <gsl/gsl_sf_gamma.h>
#include "gk15.h"
#include <errno.h>
#include <math.h>
#include <error.h>


/* UDF use real precision */
#if RP_DOUBLE == 0
#error "UDF requires real precision mode"
#endif

typedef struct ParametryFce
{
    real eps;
    real alpha;
    real d_1;

} ParametryFce;

real IntegraceF(real, void*);


#ifdef ERRNO_CHECKING

double exps(double a, const char* file, int line)
{
    double res = exp(a);

    if(errno ==  EDOM || errno == ERANGE)
    {
        error(EXIT_FAILURE, errno, "\nFile: %s:%i\nError in exp() returned: %e exponent: %e\nError", file, line, res, a);
    }

    return res;
}

double sqrts(double a, const char* file, int line)
{
    double res = sqrt(a);

    if(errno ==  EDOM || errno == ERANGE)
    {
        error(EXIT_FAILURE, errno, "\nFile: %s:%i\nError in sqrt() returned: %e base: %e\nError", file, line, res, a);
    }

    return res;
}

double pows(double a, double b, const char* file, int line)
{
    double res = pow(a, b);

    if(errno ==  EDOM || errno == ERANGE)
    {
        error(EXIT_FAILURE, errno, "\nFile: %s:%i\nError in pow() returned: %e base: %e exponent: %e\nError", file, line, res, a, b);
    }

    return res;
}

double gammas(double a, double b, const char* file, int line)
{
    double res = gsl_sf_gamma_inc(a, b);

    if(errno ==  EDOM || errno == ERANGE)
    {
        if(isfinite(res))
            errno = 0; /* GSL handles its errrors */
        else
            error(EXIT_FAILURE, errno, "\nFile: %s:%i\nError in gsl_sf_gamma_inc() returned: %e arg1: %e arg2: %e\nError", file, line, res, a, b);
    }

    return res;
}

#define exp(a) exps(a, __FILE__, __LINE__)
#define sqrt(a) sqrts(a, __FILE__, __LINE__)
#define pow(a, b) pows(a, b, __FILE__, __LINE__)
#define gsl_sf_gamma_inc(a, b) gammas(a, b, __FILE__, __LINE__)

#define CHECK_ERRNO \
if(errno ==  EDOM || errno == ERANGE) \
{ \
   error(EXIT_FAILURE, errno, "\nFile: %s:%i\nError", __FILE__, __LINE__); \
}
#else
#define CHECK_ERRNO
#endif

#endif
