#include "Integrator.hpp"

void Integrator::Runge23(std::vector<ExpData>* dataSet, double* modData, double A, double E, double NS, double yinf, int component)
{
    double k1, k2;
    double step, tn, yn;

    modData[0] = ExpData::EqWeight[(int)ExpData::WoodType()][OptimizationEngine::GetNumberOfParameterGroups()-1][component]*(*dataSet)[0].MassFrac(); //initial condition

    for(unsigned int i = 1; i < dataSet->size(); i++)
    {
        tn = (*dataSet)[i-1].Time();
        yn = modData[i-1];

        step = (*dataSet)[i].Time() - tn;

        k1 = ModEq::Eq(tn, yn, A, E, NS, yinf);
        k2 = ModEq::Eq(tn + (2./3.)*step, yn + (2./3.)*step*k1, A, E, NS, yinf);

        modData[i] = yn + step*((1./4.)*k1 + (3./4.)*k2);

    }
}
