#ifndef LUO_MB_H_INCLUDED
#define LUO_MB_H_INCLUDED

#include "udf.h"
#include "sg_pb.h"
#include "sg_mphase.h"
#include <errno.h>
#include <error.h>
#include "gk15.h"

/* GSL uses double precision */
#if RP_DOUBLE == 0
#error "UDF requires double precision mode"
#endif

typedef struct ParametryFce
{
    real lambda;
} ParametryFce;

real DstarInt(real, void*);

#ifdef ERRNO_CHECKING

double exps(double a, const char* file, int line)
{
    double res = exp(a);

    if(errno != 0)
    {
        error(EXIT_FAILURE, errno, "\nFile: %s:%i\nError in exp() exponent: %f\nError", file, line, a);
    }

    return res;
}

double sqrts(double a, const char* file, int line)
{
    double res = sqrt(a);

    if(errno != 0)
    {
        error(EXIT_FAILURE, errno, "\nFile: %s:%i\nError in sqrt() base: %f\nError", file, line, a);
    }

    return res;
}

double pows(double a, double b, const char* file, int line)
{
    double res = pow(a, b);

    if(errno != 0)
    {
        error(EXIT_FAILURE, errno, "\nFile: %s:%i\nError in pow() base: %f exponent: %f\nError", file, line, a, b);
    }

    return res;
}

#define exp(a) exps(a, __FILE__, __LINE__)
#define sqrt(a) sqrts(a, __FILE__, __LINE__)
#define pow(a, b) pows(a, b, __FILE__, __LINE__)

#define CHECK_ERRNO \
if(errno != 0) \
{ \
   error(EXIT_FAILURE, errno, "\nFile: %s:%i\nError", __FILE__, __LINE__); \
}
#else
#define CHECK_ERRNO
#endif

#endif
