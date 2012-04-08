#ifndef OPTIMIZATIONENGINE_HPP_INCLUDED
#define OPTIMIZATIONENGINE_HPP_INCLUDED

#include <eo>
#include <es.h>
#include <vector>

#include "ExpData.hpp"
#include "Integrator.hpp"

typedef eoReal<eoMinimizingFitness> Indi;



class OptimizationEngine
{
private:
    static const int SEED = 1337;
    static const int POP_SIZE = 50; // Size of population
    static const int MAX_GEN = 50; // Maximum number of generation before STOP

    //initial guesses
    static const double A_initial;
    static const double E_initial;
    static const double NS_initial;
    static const double yinf_initial;


    static double FitnessFce(const std::vector<double>& pars);
    static std::vector<ExpData>* DataSet;
    static double* ModData;
    static void InitPop(eoPop<Indi>& pop, eoEvalFuncPtr<Indi, double, const std::vector<double>& > eval);

public:
    static void Run(std::vector<ExpData>* dataSet);

};



#endif // OPTIMIZATIONENGINE_HPP_INCLUDED
