#ifndef LUO_H_INCLUDED
#define LUO_H_INCLUDED

#include "udf.h"
#include "metric.h"
#include "sg_pb.h"
#include "sg_mphase.h"
#include "gk15.h"
#include "gamma_inc.h"
#include <math.h>


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

    if(!isfinite(res))
    {
        fprintf(stderr, "\nFile: %s:%i\nError in exp() returned: %e exponent: %e\n", file, line, res, a);
        abort();
    }

    return res;
}

double sqrts(double a, const char* file, int line)
{
    double res = sqrt(a);

    if(!isfinite(res))
    {
        fprintf(stderr, "\nFile: %s:%i\nError in sqrt() returned: %e base: %e\n", file, line, res, a);
        abort();
    }

    return res;
}

double pows(double a, double b, const char* file, int line)
{
    double res = pow(a, b);

    if(!isfinite(res))
    {
        fprintf(stderr, "\nFile: %s:%i\nError in pow() returned: %e base: %e exponent: %e\n", file, line, res, a, b);
        abort();
    }

    return res;
}

double gammas(double a, double b, int index, const char* file, int line)
{
    double res = gamma_inc(a, b, index);

    if(!isfinite(res))
    {

        fprintf(stderr, "\nFile: %s:%i\nError in gsl_sf_gamma_inc() returned: %e arg1: %e arg2: %e\n", file, line, res, a, b);
        abort();
    }

    return res;
}

#define exp(a) exps(a, __FILE__, __LINE__)
#define sqrt(a) sqrts(a, __FILE__, __LINE__)
#define pow(a, b) pows(a, b, __FILE__, __LINE__)
#define gsl_sf_gamma_inc(a, b, c) gammas(a, b, c, __FILE__, __LINE__)

#define CHECK_ERRNO(result) \
if(!isfinite(result)) { \
   fprintf(stderr, "\nDetected: %e in %s:%i\n",  result, __FILE__, __LINE__); \
   abort(); } \
   (void)0
#else
#define CHECK_ERRNO (void)0
#endif

#endif
