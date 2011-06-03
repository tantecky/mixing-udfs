#ifndef LUO_MB_H_INCLUDED
#define LUO_MB_H_INCLUDED

#include "udf.h"
#include "sg_pb.h"
#include "sg_mphase.h"
#include <errno.h>
#include <error.h>
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
