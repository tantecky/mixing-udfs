#include "LuoLite.h"

#define INTEGRACNI_LIMIT 1000
#define INTEGRACNI_CHYBA_ABS 1e-7
#define INTEGRACNI_CHYBA_REL 0.0
/*fyzikalni konstanty*/
#define SIGMA 0.0728
#define RHO_L 998.2
#define C1 1.
#define RHO_G 1.225
#define MJU_L 0.001003 /*dynamicka viskozita*/


static gsl_integration_workspace* workspace = NULL;

DEFINE_EXECUTE_AT_EXIT(uvolneni_integratoru)
{
#if !RP_HOST
    if(workspace != NULL)
    {
        gsl_integration_workspace_free(workspace);

        workspace = NULL;
    }


#endif
}

DEFINE_EXECUTE_ON_LOADING(inicializace_integratoru, libname)
{
#if !RP_HOST
    if(workspace == NULL)
        workspace = gsl_integration_workspace_alloc(INTEGRACNI_LIMIT);

    if(workspace != NULL)
    {
        Message("Integracni knihovna alokovana\n");
    }
    else
    {
        Message("Nepodarilo se alokovat integracni knihovnu\n");
        abort();
    }
#endif
}

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

real IntegraceKsi(real ksi, void* parametry)
{
    ParametryFce* pars = (ParametryFce*)parametry;

    return 0.9238*pow(pars->eps, 1./3.)*pow(pars->d_1, -2./3.)*(1 - pars->alpha)\
    *pow((1 + ksi), 2)*pow(ksi, -11./3.)\
           *exp(-12*pars->cf*SIGMA*pow(pars->eps, -2./3.)*pow(pars->d_1, -5./3.)*pow(ksi, -11./3.)*(1./RHO_L));

}

DEFINE_PB_BREAK_UP_RATE_PDF(break_up_pdf_par, cell, thread, d_1, d_2)
{
    real eps = C_D(cell, THREAD_SUPER_THREAD(thread));
    real alpha = C_VOF(cell, thread);
    real f = pow(d_2, 2.0)*pow(d_1, -2.0);

    ParametryFce pars = {eps, alpha, d_1, pow(f, 2./3.) + pow(1 - f, 2./3.) - 1.}; /*0 je cfd ktere se bude menit behem integrace*/

    gsl_function fce;
    fce.function = &IntegraceKsi;
    fce.params = &pars;

    real result, error;

    real b = 1;
    real eta = pow(pow(MJU_L/RHO_L, 3)/eps, 1./4.);
    real ksimin = 11.4*eta/d_1;


    int status = gsl_integration_qags(&fce, ksimin, b, INTEGRACNI_CHYBA_ABS, INTEGRACNI_CHYBA_REL, INTEGRACNI_LIMIT, workspace, &result, &error);

    if(status != GSL_SUCCESS)
    {
        Message("UDF integrace se nezdarila\nChyba: %s\n", gsl_strerror(status));
        abort();
    }

    return result;

}
