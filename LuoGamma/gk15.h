#ifndef GK_H_INCLUDED
#define GK_H_INCLUDED

#include <errno.h>
#include <math.h>
#include <error.h>
#include <stdlib.h>

#define CHECK_ERRNO \
if(errno ==  EDOM || errno == ERANGE) \
{ \
   error(EXIT_FAILURE, errno, "\nFile: %s:%i\nError", __FILE__, __LINE__); \
}

static const double nodes[7] =
{
    0.991455371120813,
    0.949107912342759,
    0.864864423359769,
    0.741531185599394,
    0.586087235467691,
    0.405845151377397,
    0.207784955007898
};

static const double weights[8] =
{
    0.022935322010529,
    0.063092092629979,
    0.104790010322250,
    0.140653259715525,
    0.169004726639267,
    0.190350578064785,
    0.204432940075298,
    0.209482141084728
};

double gk15( double (*)(double, void*), double, double, void*);

#endif
