//Integrator: explicit second-order Runge-Kutta method with two stages
#include "Integrator.hpp"

void Integrator::Runge23(std::vector<ExpData>* dataSet, double* modData, double preFac, double actEner, double reacOrd, double massFracInf, int component)
{
    double increment1, increment2;
    //double increment1, increment2, increment3, increment4;
    double step, time_N, massFrac_N;

    modData[0] = ExpData::EqWeight[(int)ExpData::WoodType()][OptimizationEngine::GetNumberOfParameterGroups()-1][component]*(*dataSet)[0].MassFrac(); //initial condition

    for(unsigned int i = 1; i < dataSet->size(); i++)
    {
        time_N = (*dataSet)[i-1].Time();
        massFrac_N = modData[i-1];

        step = (*dataSet)[i].Time() - time_N;

        increment1 = ModEq::Equation(time_N, massFrac_N, preFac, actEner, reacOrd, massFracInf);
        increment2 = ModEq::Equation(time_N + (2./3.)*step, massFrac_N + (2./3.)*step*increment1, preFac, actEner, reacOrd, massFracInf);

        modData[i] = massFrac_N + step*((1./4.)*increment1 + (3./4.)*increment2);


        /*increment1 = step*ModEq::Equation(time_N, modData[i-1], preFac, actEner, reacOrd, massFracInf);
        increment2 = step*ModEq::Equation(time_N + (1./2.)*step, modData[i-1] + (1./2.)*increment1, preFac, actEner, reacOrd, massFracInf);
        increment3 = step*ModEq::Equation(time_N + (1./2.)*step, modData[i-1] + (1./2.)*increment2, preFac, actEner, reacOrd, massFracInf);
        increment4 = step*ModEq::Equation(time_N + step, modData[i-1] + increment3, preFac, actEner, reacOrd, massFracInf);

        modData[i] = massFrac_N + 1./6.*(increment1 + 2.0*increment2 + 2.0*increment3 + increment4);*/

    }
}
