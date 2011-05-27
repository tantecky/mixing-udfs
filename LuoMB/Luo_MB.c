#include "Luo_MB.h"

#define INTEGRACNI_CHYBA_ABS 0
#define INTEGRACNI_CHYBA_REL 1e-2
/*fyzikalni konstanty*/
#define SIGMA 0.0728
#define C 0.3
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

DEFINE_PB_BREAK_UP_RATE_FREQ(break_up_freq_martinez_bazan, cell, thread, d_1)
{
    real eps = C_D(cell, THREAD_SUPER_THREAD(thread));

    return 0.25*sqrt(8.2*pow(eps*d_1,2./3.)-12.*SIGMA/RHO_L/d_1)/d_1;
}

DEFINE_PB_BREAK_UP_RATE_PDF(break_up_pdf_par, cell, thread, d_1, d_2)
{
    real eps = C_D(cell, THREAD_SUPER_THREAD(thread));
    real Dstar = d_2/d_1;
    real dc = pow(12*SIGMA/8.2/RHO_L,3./5.)/pow(eps,2./5.);
    real lambda = dc/d_1;
    real dmin = pow(12*SIGMA/8.2/RHO_L/d_1,3./2.)/eps;
    real dmax = d_1*pow(1-pow(dmin/d_1,3.),1./3.);
    real Dstarmin = dmin/d_1;
    real Dstarmax = dmax/d_1;
    real cit = (pow(Dstar,2./3.)-pow(lambda,5./3.))*(pow(1-pow(Dstar,3.),2./9.)-pow(lambda,5./3.));

    ParametryFce pars = {lambda};

    gsl_function fce;
    fce.function = &DstarInt;
    fce.params = &pars;
    real result, error;
    size_t order;

    int status = gsl_integration_qng(&fce, Dstarmin, Dstarmax, INTEGRACNI_CHYBA_ABS, INTEGRACNI_CHYBA_REL, &result, &error, &order);

    if(status != GSL_SUCCESS)
    {
        Message("UDF integrace se nezdarila\nChyba: %s\n", gsl_strerror(status));
        abort();
    }

    return cit/result;
}

real DstarInt(real Dstar, void* parametry)
{
    ParametryFce* pars = (ParametryFce*)parametry;

    return (pow(Dstar,2./3.)-pow(pars->lambda,5./3.))*(pow(1-pow(Dstar,3.),2./9.)-pow(pars->lambda,5./3.));
}

DEFINE_EXCHANGE_PROPERTY(schiller_modified,cell,mix_thread,s_col,f_col)
{
    Thread *thread_l, *thread_g;
    real x_vel_l, y_vel_l, z_vel_l, x_vel_g, y_vel_g, z_vel_g, slip_x, slip_y, slip_z, rho_l, rho_g, alpha_l, alpha_g, mu_l, mu_t_l, d_g, abs_v, rey_p, rey_p_mod, c_D, f, tau_g, K_gl;


    /* find the threads for the liquid (primary) */
    /* and gas (secondary phases)                */

    thread_l = THREAD_SUB_THREAD(mix_thread, s_col); /*liquid phase*/
    thread_g = THREAD_SUB_THREAD(mix_thread, f_col); /*gas phase*/

    /*find phase velocities and properties*/

    x_vel_l = C_U(cell, thread_l);
    y_vel_l = C_V(cell, thread_l);
    z_vel_l = C_W(cell, thread_l);

    x_vel_g = C_U(cell, thread_g);
    y_vel_g = C_V(cell, thread_g);
    z_vel_g = C_W(cell, thread_g);

    slip_x = x_vel_l - x_vel_g;
    slip_y = y_vel_l - y_vel_g;
    slip_z = z_vel_l - z_vel_g;

    rho_l = C_R(cell, thread_l);
    rho_g = C_R(cell, thread_g);

    alpha_l = C_VOF(cell, thread_l);
    alpha_g = C_VOF(cell, thread_g);

    mu_l = C_MU_L(cell, thread_l);
    mu_t_l = C_MU_T(cell, thread_l);

    d_g = C_PHASE_DIAMETER(cell,thread_g);

    /*compute slip*/
    abs_v = sqrt(slip_x*slip_x + slip_y*slip_y + slip_z*slip_z);

    /*compute Reynold's number*/
    rey_p = rho_l*abs_v*d_g/mu_l + 1e-5;
    /*modified Reynold's number - used only in the drag coefficient correlation!!!*/
    rey_p_mod = rho_l*abs_v*d_g/(mu_l+C*mu_t_l) + 1e-5;

    /*compute drag and return drag exchange coeff K_gl*/
    if(rey_p_mod<=1000.)
        c_D = 24.*(1.+0.15*pow(rey_p_mod,0.687))/rey_p_mod;
    else
        c_D = 0.44;

    f = c_D*rey_p/24.;

    tau_g = rho_g*pow(d_g,2.)/(18.*mu_l);

    K_gl = alpha_l*alpha_g*rho_g*f/tau_g;

    return K_gl;
}
