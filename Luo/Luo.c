/*if ERRNO_CHECKING is defined than CHECK_ERRNO takes effect also pow, sqrt, pow are using ERRNO checking*/
#define ERRNO_CHECKING
#include "Luo.h"

/*fyzikalni konstanty*/
#define SIGMA 0.0728
#define RHO_L 998.2
#define C1 1.
#define RHO_G 1.225
#define MJU_L 0.001003 /*dynamicka viskozita*/

DEFINE_PB_COALESCENCE_RATE(aggregation_kernel_luo,cell,thread,d_1,d_2)
{
    real eps = C_D(cell, THREAD_SUPER_THREAD(thread));

    real u_1 = 1.43*pow(eps*d_1,1./3.);
    real u_2 = 1.43*pow(eps*d_2,1./3.);
    real u_12 = pow(pow(u_1,2.)+pow(u_2,2.),1./2.);
    real x_12 = d_1/d_2;
    real omega = M_PI/4.*pow((d_1+d_2),2.)*u_12;
    real We = RHO_L*d_1*u_12*u_12*(1./SIGMA);
    real Pag = exp(-C1*sqrt(0.75*(1+x_12*x_12)*(1+x_12*x_12*x_12))*pow(RHO_G/RHO_L + 0.5, -1./2.)*pow(1+x_12, -3.)*sqrt(We));

    real res = omega*Pag;

    CHECK_ERRNO

    return res;
}

DEFINE_PB_BREAK_UP_RATE_FREQ(break_up_freq_luo, cell, thread, d_1)
{
    real eps = C_D(cell, THREAD_SUPER_THREAD(thread));
    real alpha = C_VOF(cell, thread);

    ParametryFce pars = {eps, alpha, d_1, 0.}; /*0 je cfd ktere se bude menit behem integrace*/

    real a = 0;
    real b = 1;

    real result = gk15(&IntegraceF, a, b, &pars);

    CHECK_ERRNO

    return 0.5*result;
}

real IntegraceKsi(real ksi, void* parametry)
{
    ParametryFce* pars = (ParametryFce*)parametry;

    return 0.9238*pow(pars->eps, 1./3.)*pow(pars->d_1, -2./3.)*(1 - pars->alpha)*pow((1 + ksi), 2)*pow(ksi, -11./3.)\
           *exp(-12*(pow(pars->f, 2./3.) + pow(1 - pars->f, 2./3.) - 1)*SIGMA*pow(pars->eps, -2./3.)*pow(pars->d_1, -5./3.)*pow(ksi, -11./3.)*(1./RHO_L));

}

real IntegraceF(real f, void* parametry)
{
    ParametryFce* pars = (ParametryFce*)parametry;
    pars->f = f;

    real b = 1;
    real eta = pow(pow(MJU_L/RHO_L, 3)/pars->eps, 1./4.);
    real ksimin = 11.4*eta/pars->d_1;

    CHECK_ERRNO

    real result = gk15(&IntegraceKsi, ksimin, b, pars);

    CHECK_ERRNO

    return result;
}

DEFINE_PB_BREAK_UP_RATE_PDF(break_up_pdf_par, cell, thread, d_1, d_2)
{
    real eps = C_D(cell, THREAD_SUPER_THREAD(thread));
    real alpha = C_VOF(cell, thread);

    ParametryFce pars = {eps, alpha, d_1, 0.}; /*0 je cfd ktere se bude menit behem integrace*/


    real a = 0;
    real b = 1;

    real result = gk15(&IntegraceF, a, b,  &pars);

    CHECK_ERRNO

    /*jmenovatel*/

    real betaDen = 0.5*result;

    pars.f = pow(d_2, 3.)*pow(d_1, -3.);

    b = 1;
    real eta = pow(pow(MJU_L/RHO_L, 3)/eps, 1./4.);
    real ksimin = 11.4*eta/d_1;

    CHECK_ERRNO

    result = gk15(&IntegraceKsi, ksimin, b, &pars);

    real res = result/betaDen;

    CHECK_ERRNO

    return res;
}
