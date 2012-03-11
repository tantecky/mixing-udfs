#ifndef _PRINCE_H_
#define _PRINCE_H_

#include "udf.h"
#include "udf.h"
#include "sg_pb.h"
#include "sg_mphase.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

/* GSL je cela double */
#if RP_DOUBLE == 0
#error "UDF vyzaduje double precision mode"
#endif

#ifdef DEBUG_VYPIS
#if PARALLEL
#error "DEBUG_VYPIS vyzaduje serial resic"
#endif
#endif

typedef struct ParametryFce
{
    real eps;
    real alpha;
    real d_1;

} ParametryFce;

real IntegrovanaFce(real, void *);

#endif
