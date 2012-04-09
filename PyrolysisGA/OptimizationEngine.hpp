#ifndef OPTIMIZATIONENGINE_HPP_INCLUDED
#define OPTIMIZATIONENGINE_HPP_INCLUDED

#include <eo>
#include <es.h>
#include <vector>

#include "ExpData.hpp"
#include "Integrator.hpp"

typedef eoReal<eoMaximizingFitness> Indi;



class OptimizationEngine
{
private:
    static const int SEED;
    static const int POP_SIZE; // Size of population
    static const int MAX_GEN; // Maximum number of generation before STOP

    static const double ALFA;     //BLX coefficient
    static const double HYPER_CUBE_RATE;     // relative weight for hypercube Xover
    static const double SEGMENT_RATE;  // relative weight for segment Xover
    static const double EPSILON;	// range for real uniform mutation
    static const double SIGMA;	    	// std dev. for normal mutation
    static const double UNIFORM_MUT_RATE;  // relative weight for uniform mutation
    static const double DET_MUT_RATE;      // relative weight for det-uniform mutation
    static const double NORMAL_MUT_RATE;  // relative weight for normal mutation
    static const double P_CROSS;	// Crossover probability
    static const double P_MUT;	// mutation probability

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

