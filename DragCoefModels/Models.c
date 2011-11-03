/*
Several drag models for suspension
    - Pinelli
    - Brucato
    - Khopkar

Author: Tomáš Antecký (21. 10. 2011)
Develop for GCC according to gnu90 standard
Based on André Bakker code
*/

#include <math.h>
#include <udf.h>
#include <sg_mphase.h>
#include <sg.h>

#define DIAMETER 1.02e-3 /*diameter of solid particle*/
#define RHO_S 1136.02 /*density of solid phase*/
#define RHO_L 1011.44 /*density of liquid phase*/
#define MU_L 5e-3 /*dynamic viscosity of liquid phase*/
#define NU_L (MU_L/RHO_L) /*kinematic viscosity of liquid phase*/

/*checking if exchange coefficient is finite number*/
#define DEBUG_COEF

#ifdef DEBUG_COEF
#define CHECK_COEF(result) \
if(!isfinite(result)) { \
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

    /*drag force*/
    real fdrag = cd*reyp/24.;

    /*volumetric fraction of liquid phase*/
    real vol_l = C_VOF(cell, thread_l);

    real k_s_l = (1. - vol_l)*vol_l*RHO_S*fdrag/taup;

    CHECK_COEF(k_s_l);

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

    /*drag force*/
    real fdrag = cd*reyp/24.;

    /*volumetric fraction of liquid phase*/
    real vol_l = C_VOF(cell, thread_l);

    real k_s_l = (1. - vol_l)*vol_l*RHO_S*fdrag/taup;

    CHECK_COEF(k_s_l);

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

    /*drag force*/
    real fdrag = cd*reyp/24.;

    /*volumetric fraction of liquid phase*/
    real vol_l = C_VOF(cell, thread_l);

    real k_s_l = (1. - vol_l)*vol_l*RHO_S*fdrag/taup;

    if(isnan(k_s_l))
        return 0.0;

    CHECK_COEF(k_s_l);

    return k_s_l;

}
