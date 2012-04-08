#ifndef INTEGRATOR_HPP_INCLUDED
#define INTEGRATOR_HPP_INCLUDED

#include <vector>
#include "ModEq.hpp"
#include "ExpData.hpp"

class Integrator
{

public:
    static void Runge23(std::vector<ExpData>* dataSet, double* modData, double A, double E, double NS, double yinf);

};

#endif // INTEGRATOR_HPP_INCLUDED
