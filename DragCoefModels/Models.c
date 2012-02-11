/*
Several drag coefficient models for solid-liquid interphase drag force
    - Schiller-Nauman
    - Pinelli
    - Brucato
    - Khopkar

Author: Tomáš Antecký (21. 10. 2011)
Develop for GCC according to gnu90 standard
*/

#include <math.h>
#include <udf.h>
#include <sg_mphase.h>
#include <sg.h>

#define DIAMETER (1.02e-3) /*diameter of the solid particle*/
#define RHO_L (1011.44) /*density of the liquid phase*/
#define MU_L (5e-3) /*dynamic viscosity of the liquid phase*/
#define NU_L (MU_L/RHO_L) /*kinematic viscosity of the liquid phase*/

/*checking if exchange coefficient is a finite number*/
#define DEBUG_COEF

/*if is defined the drag coeffcient is stored in user defined memory*/
#define USE_UDM

#ifdef USE_UDM
#define UDM_NUM (0) /*the id of the drag coefficient in UDM*/
#endif

#ifdef DEBUG_COEF
#define CHECK_COEF(result) \
if(!isfinite(result) || result < 0.0) { \
   fprintf(stderr, "\nDetected: %e in %s:%i\n", result, __FUNCTION__, __LINE__); \
   abort(); } \
   (void)0
#else
#define CHECK_COEF(result) (void)0
#endif

DEFINE_EXECUTE_ON_LOADING(on_load, libname)
{
    Message("\nBuilded: %s %s\n", __DATE__, __TIME__);
}

DEFINE_EXCHANGE_PROPERTY(SchillerNauman_CD, cell, mix_thread, s_col, f_col)
{
    /*cd0 is Schiller-Nauman's drag coefficinet*/
#include "SchillerNauman.h"

    /*volumetric fraction of the solid phase*/
    real vol_s = C_VOF(cell, thread_s);

    real k_s_l = (3./4.)*vol_s/DIAMETER*cd0*RHO_L*slip;

    if(isnan(k_s_l))
        return 0.0;

    CHECK_COEF(k_s_l);

#ifdef USE_UDM
    C_UDMI(cell, mix_thread, UDM_NUM) = cd0;
#endif

    return k_s_l;

}

DEFINE_EXCHANGE_PROPERTY(Pinelli_CD, cell, mix_thread, s_col, f_col)
{
    /*cd0 is Schiller-Nauman's drag coefficinet*/
#include "SchillerNauman.h"

    /*---Pinelli Correction---*/

    /*energy dissipation rate*/
    real eps = C_D(cell,thread_l);

    /*Kolmogoroff microscale*/
    real kolscale = pow((NU_L*NU_L*NU_L)/eps, 0.25);

    real pinelli = 0.4*tanh(16.*kolscale/DIAMETER - 1.) + 0.6;

    real cd = cd0/pow(pinelli, 2.0);

    /*if (reyp<0.001)
        reyp=0.001;*/

    /*volumetric fraction of the solid phase*/
    real vol_s = C_VOF(cell, thread_s);

    real k_s_l = (3./4.)*vol_s/DIAMETER*cd*RHO_L*slip;

    if(isnan(k_s_l))
        return 0.0;

    CHECK_COEF(k_s_l);

#ifdef USE_UDM
    C_UDMI(cell, mix_thread, UDM_NUM) = cd;
#endif

    return k_s_l;

}

DEFINE_EXCHANGE_PROPERTY(Brucato_CD, cell, mix_thread, s_col, f_col)
{
    /*cd0 is Schiller-Nauman's drag coefficinet*/
#include "SchillerNauman.h"

    /*---Brucato Correction---*/

    /*energy dissipation rate*/
    real eps = C_D(cell,thread_l);

    /*Kolmogoroff microscale*/
    real kolscale = pow((NU_L*NU_L*NU_L)/eps, 0.25);

    real cd = cd0*(1 + 8.76e-4*pow(DIAMETER/kolscale, 3.0));

    /*if (reyp<0.001)
        reyp=0.001;*/

    /*volumetric fraction of the solid phase*/
    real vol_s = C_VOF(cell, thread_s);

    real k_s_l = (3./4.)*vol_s/DIAMETER*cd*RHO_L*slip;

    if(isnan(k_s_l))
        return 0.0;

    CHECK_COEF(k_s_l);

#ifdef USE_UDM
    C_UDMI(cell, mix_thread, UDM_NUM) = cd;
#endif

    return k_s_l;

}

DEFINE_EXCHANGE_PROPERTY(Khopkar_CD, cell, mix_thread, s_col, f_col)
{
    /*cd0 is Schiller-Nauman's drag coefficinet*/
#include "SchillerNauman.h"

    /*---Khopkar Correction---*/

    /*energy dissipation rate*/
    real eps = C_D(cell,thread_l);

    /*Kolmogoroff microscale*/
    real kolscale = pow((NU_L*NU_L*NU_L)/eps, 0.25);

    real cd = cd0*(1 + 8.76e-5*pow(DIAMETER/kolscale, 3.0));

    /*if (reyp<0.001)
        reyp=0.001;*/

    /*volumetric fraction of the solid phase*/
    real vol_s = C_VOF(cell, thread_s);

    real k_s_l = (3./4.)*vol_s/DIAMETER*cd*RHO_L*slip;

    if(isnan(k_s_l))
        return 0.0;

    CHECK_COEF(k_s_l);

#ifdef USE_UDM
    C_UDMI(cell, mix_thread, UDM_NUM) = cd;
#endif

    return k_s_l;
}

DEFINE_EXCHANGE_PROPERTY(Debug_CD, cell, mix_thread, s_col, f_col)
{
    printf("s_col: %d\n", s_col);
    printf("f_col: %d", f_col);
    abort();

    return 0;

}

DEFINE_ON_DEMAND(Quality_of_suspension)
{
    Domain* d;
    Thread *t;
    cell_t c;
    d = Get_Domain(1);
    /* Get the domain using ANSYS FLUENT utility */
    /* Loop over all cell threads in the domain */

    real maxFrac = -1;
    real frac;

    thread_loop_c(t,d)
    {
        /* Compute max, min, volume-averaged temperature */
        /* Loop over all cells */
        begin_c_loop(c,t)
        {
            frac = C_VOF(c, THREAD_SUB_THREAD(t, 1));

            if(frac > maxFrac)
                maxFrac = frac;
        }
        end_c_loop(c,t)

    }
    Message0("maxFrac: %e", maxFrac);
}
