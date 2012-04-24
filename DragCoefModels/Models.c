/*
Several drag coefficient models for solid-liquid interphase drag force
    - Schiller-Nauman
    - Pinelli
    - Brucato
    - Khopkar

Author: Tomáš Antecký (21. 10. 2011)
Developed for GCC and MS compiler
Tested on Fluent 12.1.4
*/

#include <math.h>
#include <udf.h>
#include <sg_mphase.h>
#include <sg.h>

#ifdef _MSC_VER /*for MS compiler*/
#include <float.h>
#define isfinite(x) _finite(x)
#define isnan(x) _isnan(x)
#endif

/*======SETTINGS STARTS HERE======*/
/*physical properties*/
#define DIAMETER (1.02e-3) /*diameter of the solid particle*/
#define RHO_L (999.5) /*density of the liquid phase*/
#define MU_L (1.138e-3) /*dynamic viscosity of the liquid phase*/
#define NU_L (MU_L/RHO_L) /*kinematic viscosity of the liquid phase*/

#define SOLID_PHASE_ID (1) /*id of the solid phase*/

/*checking if computed exchange coefficient is a finite number*/
#define DEBUG_COEF

/*if is defined the drag coeffcient is stored in user defined memory*/
/*don't forget to allocate memory for it under:
  Define -> User-Defined -> Memory*/
#define USE_UDM

#ifdef USE_UDM
#define CD_UDM_NUM (0) /*the id of the drag coefficient in UDM*/
#endif
/*======SETTINGS ENDS HERE======*/

#ifdef DEBUG_COEF
#define CHECK_COEF(result) \
    if(!isfinite(result) || result < 0.0) { \
       fprintf(stderr, "\nDetected: %e in %s:%i\n", result, __FUNCTION__, __LINE__); \
       abort(); } \
       (void)0
#else
#define CHECK_COEF(result) (void)0
#endif

/*prototypes*/
real QualityOfSuspension(void);
real MeanVolFrac(void);

DEFINE_EXECUTE_ON_LOADING(on_load, libname)
{
    Message("\nBuilded: %s %s\n", __DATE__, __TIME__);
}

DEFINE_EXCHANGE_PROPERTY(SchillerNauman_CD, cell, mix_thread, s_col, f_col)
{
    real vol_s;
    real k_s_l;

    /*cd0 is Schiller-Nauman's drag coefficinet*/
#include "SchillerNauman.h"

    /*volumetric fraction of the solid phase*/
    vol_s = C_VOF(cell, thread_s);

    k_s_l = (3./4.)*vol_s/DIAMETER*cd0*RHO_L*slip;

    if(isnan(k_s_l))
        return 0.0;

    CHECK_COEF(k_s_l);

#ifdef USE_UDM
    C_UDMI(cell, mix_thread, CD_UDM_NUM) = cd0;
#endif

    return k_s_l;

}

DEFINE_EXCHANGE_PROPERTY(Pinelli_CD, cell, mix_thread, s_col, f_col)
{
    real eps;
    real kolscale;
    real cd;
    real vol_s;
    real k_s_l;
    real pinelli;

    /*cd0 is Schiller-Nauman's drag coefficinet*/
#include "SchillerNauman.h"

    /*---Pinelli Correction---*/

    /*energy dissipation rate*/
    eps = C_D(cell,thread_l);

    /*Kolmogoroff microscale*/
    kolscale = pow((NU_L*NU_L*NU_L)/eps, 0.25);

    pinelli = 0.4*tanh(16.*kolscale/DIAMETER - 1.) + 0.6;

    cd = cd0/pow(pinelli, 2.0);

    /*volumetric fraction of the solid phase*/
    vol_s = C_VOF(cell, thread_s);

    k_s_l = (3./4.)*vol_s/DIAMETER*cd*RHO_L*slip;

    if(isnan(k_s_l))
        return 0.0;

    CHECK_COEF(k_s_l);

#ifdef USE_UDM
    C_UDMI(cell, mix_thread, CD_UDM_NUM) = cd;
#endif

    return k_s_l;

}

DEFINE_EXCHANGE_PROPERTY(Brucato_CD, cell, mix_thread, s_col, f_col)
{
    real eps;
    real kolscale;
    real cd;
    real vol_s;
    real k_s_l;

    /*cd0 is Schiller-Nauman's drag coefficinet*/
#include "SchillerNauman.h"

    /*---Brucato Correction---*/

    /*energy dissipation rate*/
    eps = C_D(cell,thread_l);

    /*Kolmogoroff microscale*/
    kolscale = pow((NU_L*NU_L*NU_L)/eps, 0.25);

    cd = cd0*(1 + 8.76e-4*pow(DIAMETER/kolscale, 3.0));

    /*volumetric fraction of the solid phase*/
    vol_s = C_VOF(cell, thread_s);

    k_s_l = (3./4.)*vol_s/DIAMETER*cd*RHO_L*slip;

    if(isnan(k_s_l))
        return 0.0;

    CHECK_COEF(k_s_l);

#ifdef USE_UDM
    C_UDMI(cell, mix_thread, CD_UDM_NUM) = cd;
#endif

    return k_s_l;

}

DEFINE_EXCHANGE_PROPERTY(Khopkar_CD, cell, mix_thread, s_col, f_col)
{
    real eps;
    real kolscale;
    real cd;
    real vol_s;
    real k_s_l;

    /*cd0 is Schiller-Nauman's drag coefficinet*/
#include "SchillerNauman.h"

    /*---Khopkar Correction---*/

    /*energy dissipation rate*/
    eps = C_D(cell,thread_l);

    /*Kolmogoroff microscale*/
    kolscale = pow((NU_L*NU_L*NU_L)/eps, 0.25);

    cd = cd0*(1 + 8.76e-5*pow(DIAMETER/kolscale, 3.0));

    /*volumetric fraction of the solid phase*/
    vol_s = C_VOF(cell, thread_s);

    k_s_l = (3./4.)*vol_s/DIAMETER*cd*RHO_L*slip;

    if(isnan(k_s_l))
        return 0.0;

    CHECK_COEF(k_s_l);

#ifdef USE_UDM
    C_UDMI(cell, mix_thread, CD_UDM_NUM) = cd;
#endif

    return k_s_l;
}

real MeanVolFrac(void)
{
    Domain* d;
    Thread *t;
    cell_t c;
    unsigned int numOfCells;
    real totalVolume;
    real sumVolFrac;
    real cellVol;

    numOfCells = 0;
    totalVolume = 0.0;
    sumVolFrac = 0.0;
    d = Get_Domain(1);

    thread_loop_c(t,d)
    {
        begin_c_loop(c,t)
        {
            numOfCells++;
            cellVol = C_VOLUME(c,t);
            totalVolume += cellVol;
            sumVolFrac += cellVol*C_VOF(c, THREAD_SUB_THREAD(t, SOLID_PHASE_ID)); /*1 - secondary phase = solid phase*/
        }
        end_c_loop(c,t)

    }

    return (sumVolFrac / totalVolume);

}

real QualityOfSuspension()
{
    Domain* d;
    Thread *t;
    cell_t c;
    real frac;
    unsigned int numOfCells;
    real avgVolFrac;
    real parcSum;

    numOfCells = 0;
    avgVolFrac = MeanVolFrac();
    parcSum = 0.0;
    d = Get_Domain(1);

    thread_loop_c(t,d)
    {
        begin_c_loop(c,t)
        {
            numOfCells++;
            frac = C_VOF(c, THREAD_SUB_THREAD(t, SOLID_PHASE_ID));

            parcSum += pow(frac/avgVolFrac - 1.0, 2.0);

        }
        end_c_loop(c,t)

    }

    return sqrt(parcSum/numOfCells);

}

/*compute and display on demand quality of suspension*/
/*Define -> User-Defined -> Execute on Demand*/
DEFINE_ON_DEMAND(QA_of_suspension)
{
    real qa;

    qa = QualityOfSuspension();

    Message0("\nQualityOfSuspension: %f\n", qa);
}
