#ifndef MODEQ_HPP_INCLUDED
#define MODEQ_HPP_INCLUDED

#include <cmath>
#include <iostream>
#include <cstdlib>
#include "ExpData.hpp"

class ModEq
{
private:
    static const double R = 8.314472;
    static const double T0 =  17 + 273.15;

public:
    static double Eq(double t, double y, double A, double E, double NS, double yinf);
};


#endif // MODEQ_HPP_INCLUDED
