#ifndef LUO_MB_H_INCLUDED
#define LUO_MB_H_INCLUDED

#include "udf.h"
#include "sg_pb.h"
#include "sg_mphase.h"
#include "gk15.h"

/* GSL use real precision */
#if RP_DOUBLE == 0
#error "UDF requires real precision mode"
#endif

typedef struct ParametryFce
{
    real lambda;
} ParametryFce;

real DstarInt(real, void*);

#endif
