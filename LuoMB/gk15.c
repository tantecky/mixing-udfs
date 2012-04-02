#include "gk15.h"

double gk15( double (*fce)(double, void*), double a, double b, void* pars)
{
    double trans1 = (b-a)/2.0;
    double trans2 = (a+b)/2.0;

    double integral = 0.0;

    int i;
    for(i = 0; i < 7; i++)
    {
        integral += weights[i]*((*fce)(trans1*nodes[i] + trans2, pars) + (*fce)(trans1*(-nodes[i]) + trans2, pars));
    }

    integral += weights[7]*(*fce)(trans2, pars);

    return trans1*integral;
}
