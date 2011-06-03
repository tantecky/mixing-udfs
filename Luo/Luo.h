#ifndef LUO_H_INCLUDED
#define LUO_H_INCLUDED

#include "udf.h"
#include "sg_pb.h"
#include "sg_mphase.h"
#include "gk15.h"
#include <errno.h>
#include <error.h>

/* GSL use real precision */
#if RP_DOUBLE == 0
#error "UDF requires real precision mode"
#endif

typedef struct ParametryFce
{
    real eps;
    real alpha;
    real d_1;
    real f;

} ParametryFce;

real IntegraceKsi(real, void*);
real IntegraceF(real, void*);

#ifdef DEBUG
#define CHECK_ERRNO \
if(errno != 0) \
{ \
   error(EXIT_FAILURE, errno, "\nFile: %s:%i\nError", __FILE__, __LINE__); \
}
#else
#define CHECK_ERRNO
#endif

#endif
