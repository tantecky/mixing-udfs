#ifndef MODEQ_HPP_INCLUDED
#define MODEQ_HPP_INCLUDED

#include <cmath>
#include <iostream>
#include <cstdlib>
#include "ExpData.hpp"

class ModEq
{
private:
    static const double UNI_GAS_CONST = 8.314462;
    static const double TEMP_INIT =  17 + 273.15;

public:
    static double Equation(double time, double massFrac, double preFac, double actEner, double reacOrd, double massFracInf);
};


#endif // MODEQ_HPP_INCLUDED
