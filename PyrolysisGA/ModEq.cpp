//Model equation: Arrhenius equation
#include "ModEq.hpp"

double ModEq::Equation(double time, double massFrac, double preFac, double actEner, double reacOrd, double massFracInf)
{
    double beta = ExpData::Beta();
    return -preFac*( pow(massFrac-massFracInf,reacOrd))*exp(-actEner/(UNI_GAS_CONST*(beta*time+TEMP_INIT)));
}
