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
#include <stdint.h>
#include <udf.h>
#include <sg_mphase.h>
#include <sg.h>

#define DIAMETER (1.02e-3) /*diameter of the solid particle*/
#define RHO_L (1011.44) /*density of the liquid phase*/
#define MU_L (5e-3) /*dynamic viscosity of the liquid phase*/
#define NU_L (MU_L/RHO_L) /*kinematic viscosity of the liquid phase*/
#define SOLID_PHASE_ID (1) /*id of the solid phase*/

/*checking if exchange coefficient is a finite number*/
#define DEBUG_COEF

/*if is defined the drag coeffcient is stored in user defined memory*/
#define USE_UDM

#ifdef USE_UDM
#define CD_UDM_NUM (0) /*the id of the drag coefficient in UDM*/
#define QA_UDM_NUM (1) /*the id of the quality of suspension in UDM*/
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

/*prototypes*/
real QualityOfSuspension(void);
real MeanVolFrac(void);

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
    C_UDMI(cell, mix_thread, CD_UDM_NUM) = cd0;
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
    C_UDMI(cell, mix_thread, CD_UDM_NUM) = cd;
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
    C_UDMI(cell, mix_thread, CD_UDM_NUM) = cd;
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
    C_UDMI(cell, mix_thread, CD_UDM_NUM) = cd;
#endif

    return k_s_l;
}

real MeanVolFrac(void)
{
    Domain* d;
    Thread *t;
    cell_t c;
    d = Get_Domain(1);

    uint_least32_t numOfCells = 0;
    real totalVolume = 0.0;
    real sumVolFrac = 0.0;
    real cellVol;

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
    d = Get_Domain(1);

    uint_least32_t numOfCells = 0;
    real avgVolFrac = MeanVolFrac();

    real frac;
    real parcSum = 0.0;

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

#ifdef USE_UDM
DEFINE_EXECUTE_AT_END(Quality_of_suspension)
{
    Domain* d;
    Thread* t;
    cell_t c;
    d = Get_Domain(1);

    real qa = QualityOfSuspension();

    Message0("\nQualityOfSuspension: %f\n", qa);

    /*thread_loop_c(t,d)
    {
        begin_c_loop(c,t)
        {
            C_UDMI(c, t, QA_UDM_NUM) = qa;
        }
        end_c_loop(c,t)

    }*/
}
#endif

DEFINE_ON_DEMAND(QA_of_suspension)
{
    real qa = QualityOfSuspension();
    real avgVolFrac = MeanVolFrac();

    Message0("\nQualityOfSuspension: %f\n", qa);
    Message0("\nMeanVolFrac: %f\n", avgVolFrac);
}
