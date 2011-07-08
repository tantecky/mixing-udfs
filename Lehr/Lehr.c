/*if ERRNO_CHECKING is defined than CHECK_ERRNO takes effect also pow, sqrt, pow are using ERRNO checking*/
#define ERRNO_CHECKING
#include "Lehr.h"

/*fyzikalni konstanty*/
#define SIGMA 0.0728
#define RHO_L 998.2
#define C1 1.
#define RHO_G 1.225
#define MJU_L 0.001003 /*dynamicka viskozita*/
#define C 0.3


DEFINE_EXECUTE_ON_LOADING(on_load, libname)
{
    Message("\nBuilded: %s %s\n", __DATE__, __TIME__);
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

    CHECK_ERRNO(Pag);

    real result = omega*Pag;

    CHECK_ERRNO(result);

    return result;
}

DEFINE_PB_BREAK_UP_RATE_FREQ(break_up_freq_luo, cell, thread, d_1)
{
    real eps = C_D(cell, THREAD_SUPER_THREAD(thread));
    real g = 0.5*pow(d_1, 5./3.)*pow(eps, 19./15.)*pow(RHO_L, 7./5.)/pow(SIGMA, 7./5.)*exp(-sqrt(2.)*pow(SIGMA, 9./5.)*pow(d_1, -3.)*pow(RHO_L, -9./5.)*pow(eps, -6./5.));

    CHECK_ERRNO(g);

    return g;
}

inline double betav(double eps, double ai, double aj)
{
    double pdf1 = 6./pow(M_PI, 3./2.)/pow(ai, 3.);
    double pdf2 = exp(-9./4.*pow(log(pow(2.,2./5.)*ai*pow(RHO_L, 3./5.)*pow(eps, 2./5.)/pow(SIGMA, 3./5.)), 2.));
    double pdf3 = 1 + erf(3./2.*log(pow(2.,1./15.)*aj*pow(RHO_L, 3./5.)*pow(eps, 2./5.)/pow(SIGMA, 3./5.)));

    CHECK_ERRNO(pdf1);
    CHECK_ERRNO(pdf2);
    CHECK_ERRNO(pdf3);

    double res = pdf1*pdf2/pdf3;

    CHECK_ERRNO(res);

    return res;
}

DEFINE_PB_BREAK_UP_RATE_PDF(break_up_pdf_par, cell, thread, d_1, d_2)
{
    real eps = C_D(cell, THREAD_SUPER_THREAD(thread));
    real ai = d_2;
    real aj = d_1;

    real f = pow(ai, 3.)/pow(aj, 3.);

    if(0. <= f && f <= 0.5)
        return betav(eps, ai, aj);
    else if(0.5 < f && f <= 1.)
        return betav(eps, pow(aj, 3.)-pow(ai, 3.), pow(aj, 3.));
    else
    {
        fprintf(stderr, "\nFile: %s:%i\nError: f is %e\n", __FILE__, __LINE__, f);
        abort();
    }
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
    if(rey_p_mod <= 1000.)
        c_D = 24.*(1.+0.15*pow(rey_p_mod,0.687))/rey_p_mod;
    else
        c_D = 0.44;

    CHECK_ERRNO(c_D);

    f = c_D*rey_p/24.;

    tau_g = rho_g*pow(d_g,2.)/(18.*mu_l);

    K_gl = alpha_l*alpha_g*rho_g*f/tau_g;

    CHECK_ERRNO(K_gl);

    return K_gl;
}

