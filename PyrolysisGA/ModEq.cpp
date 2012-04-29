#include "ModEq.hpp"

double ModEq::Eq(double t, double y, double A, double E, double NS, double yinf)
{
    double beta = ExpData::Beta();
    return -A*( pow(y-yinf,NS))*exp(-E/(R*(beta*t+T0)));
}
