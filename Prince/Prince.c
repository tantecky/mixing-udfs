#include "Prince.h"

/* NASTAVENI KONSTANT ZDE!!! */
#define INTEGRACNI_CHYBA 1e-7
#define RHO_L 998.2
#define SIGMA 0.0728
#define C1 1.0
#define H0 10.0e-4   /*zmenit pripadne na 10.0e-3*/
#define HF 10.0e-8   /*zmenit pripadne na 10.0e-6*/
#define INTEGRACNI_LIMIT 1000
/*#define DEBUG_VYPIS*/ /*vypisovat hodnoty do souboru DEBUG_VYPIS - jinak zakometovat!!!*/

static gsl_integration_workspace* workspace = NULL;

#ifdef DEBUG_VYPIS
static FILE* debug_soubor = NULL;
#endif

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
    }
#endif

#ifdef DEBUG_VYPIS
    debug_soubor = fopen("DEBUG_VYPIS", "w");

    if(debug_soubor == NULL)
    {
        Message("Nepodarilo se vytvorit soubor DEBUG_VYPIS\n");
        abort();
    }
#endif

}

DEFINE_PB_COALESCENCE_RATE(aggregation_kernel_prince,cell,thread,d_1,d_2)
{
    real agg_kernel, eps, r_1, r_2, r_12, omegaT, Pc;

    eps = C_D(cell, THREAD_SUPER_THREAD(thread));
    r_1 = d_1/2.0;
    r_2 = d_2/2.0;
    r_12 = 0.5*r_1*r_2/(r_1+r_2);
    omegaT = 0.089*3.1415*pow(eps,1./3.)*pow((d_1+d_2),2.)*pow((pow(d_1,2./3.)+pow(d_2,2./3.)),1./2.);
    Pc = exp(-pow(r_12,5./6.)*pow(RHO_L,1./2.)*pow(eps,1./3.)*log(H0/HF)/(4.*pow(SIGMA,1./2.)));

    agg_kernel = omegaT*Pc;

    return agg_kernel;
}

DEFINE_PB_BREAK_UP_RATE_FREQ(break_up_freq_prince, cell, thread, d_1)
{
    ParametryFce pars = { C_D(cell, THREAD_SUPER_THREAD(thread)), C_VOF(cell, thread), d_1 };

    gsl_function fce;
    fce.function = &IntegrovanaFce;
    fce.params = &pars;

    /* MEZE */
    real a = 0.2*d_1;
    real b = d_1;

    real result, error;

    /* popis metody:
        http://www.gnu.org/software/gsl/manual/html_node/QAGS-adaptive-integration-with-singularities.html
    */
    int status = gsl_integration_qags(&fce, a, b, INTEGRACNI_CHYBA, 0, INTEGRACNI_LIMIT, workspace, &result, &error);

    if(status != GSL_SUCCESS)
    {
        Message("UDF integrace se nezdarila\nChyba: %s\n", gsl_strerror(status));
    }

#ifdef DEBUG_VYPIS

#endif

    return C1*result;

}

real IntegrovanaFce(real ksi, void* parametry)
{
    ParametryFce* pars = (ParametryFce*)parametry;

    real f1 = pow((pars->d_1+ksi),2.)*pow((pow(pars->d_1,2./3.)+pow(ksi,2./3.)),1./2.)*pow(pars->eps,1./3.);
    real f2 = exp(-1.18/pow(ksi,2./3.)*SIGMA/(RHO_L*(pars->d_1)*pow(pars->eps,2./3.)));

    return 0.07*pow(3.1415,4.)/pow(ksi,4.)*f1*f2;
}


DEFINE_PB_BREAK_UP_RATE_PDF(break_up_pdf_par, cell, thread, d_1, d_2)
{

    return 1.0;
}
