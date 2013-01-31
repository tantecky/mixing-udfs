#ifndef INTEGRATOR_HPP_INCLUDED
#define INTEGRATOR_HPP_INCLUDED

#include <vector>
#include "ModEq.hpp"
#include "ExpData.hpp"
#include "OptimizationEngine.hpp"

class Integrator
{

public:
    static void Runge23(std::vector<ExpData>* dataSet, double* modData, double preFac, double actEner, double reacOrd, double massFracInf, int component);

};

#endif // INTEGRATOR_HPP_INCLUDED
