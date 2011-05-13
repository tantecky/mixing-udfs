#include "LuoGamma.h"

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

    return omega*Pag;
}

DEFINE_PB_BREAK_UP_RATE_FREQ(break_up_freq_luo, cell, thread, d_1)
{
    return 1.0;
}


DEFINE_PB_BREAK_UP_RATE_PDF(break_up_pdf_par, cell, thread, d_1, d_2)
{
    real eps = C_D(cell, THREAD_SUPER_THREAD(thread));
    real alpha = C_VOF(cell, thread);

    real cf = pow(d_2, 2.0)*pow(d_1, -2.0) + pow(1 - pow(d_2, 3.0)*pow(d_1, -3.0), 2./3.) - 1.;
    real beta1 = 8.0*gsl_sf_gamma(1./3.)*alpha/(5.0*M_PI);
    real k1 = 15.0*pow(M_PI, 1./3.)*sqrt(beta1)/(8. * pow(2., 2./3.)*gsl_sf_gamma(1./3.));

    real b = 12.*cf*SIGMA*pow(eps, -2./3.)*pow(d_2, -5./3.)/(beta1*MJU_L);
    real tm = b*pow(MJU_L/d_2, -11./3.);

    real g = (gsl_sf_gamma_inc(2./11., tm) - gsl_sf_gamma_inc(2./11., b));




    return g;

}
