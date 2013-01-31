#ifndef OPTIMIZATIONENGINE_HPP_INCLUDED
#define OPTIMIZATIONENGINE_HPP_INCLUDED

#include <eo>
#include <es.h>
#include <vector>
#include <sstream>
#include <cstring>
#include <fstream>
#include <iostream>
#include <utils/eoTimedMonitor.h>

#include "ExpData.hpp"
#include "Integrator.hpp"

typedef eoReal<eoMinimizingFitness> Indi;



class OptimizationEngine
{
private:
    static const int SEED;
    static const int POP_SIZE;                  // Size of population
    static const int MAX_GEN;                   // Maximum number of generation before STOP
    static const unsigned int PRINT_EVERY_SEC;  // print info to console every specified second

    static const double ALFA;                   // BLX coefficient
    static const double HYPER_CUBE_RATE;        // relative weight for hypercube Xover
    static const double SEGMENT_RATE;           // relative weight for segment Xover
    static const double EPSILON;	            // range for real uniform mutation
    static const double SIGMA;	    	        // std dev. for normal mutation
    static const double UNIFORM_MUT_RATE;       // relative weight for uniform mutation
    static const double DET_MUT_RATE;           // relative weight for det-uniform mutation
    static const double NORMAL_MUT_RATE;        // relative weight for normal mutation
    static const double P_CROSS;	            // Crossover probability
    static const double P_MUT;	                // mutation probability

    //initial guesses
    static const double PRE_FAC_INITIAL;
    static const double ACT_ENER_INITIAL;
    static const double REAC_ORD_INITIAL;
    static const double MASS_FRAC_INF_INITIAL;

    const static double PRE_FAC_MIN[3][3];
    const static double PRE_FAC_MAX[3][3];
    const static double ACT_ENER_MIN[3][3];
    const static double ACT_ENER_MAX[3][3];
    const static double REAC_ORD_MIN[3][3];
    const static double REAC_ORD_MAX[3][3];
    const static double MASS_FRAC_INF_MIN[3][3];
    const static double MASS_FRAC_INF_MAX[3][3];

    static std::vector<ExpData>* DataSet;
    static double* ModData;
    static double* IntHeap;
    static int NumberOfParameterGroups;

    static double FitnessFce(const std::vector<double>& pars);
    static void InitPop(eoPop<Indi>& pop, eoEvalFuncPtr<Indi, double, const std::vector<double>& > eval);
    static void SaveResults(double fitness, double* A, double* E, double* NS, double* yinf);
    static double Logdist(eoRng& gen, int maxExp, bool alsoNegative);
    inline static void ClearModData();
    static void ComputeModData(double A, double E, double NS, double yinf, int component);
    static std::vector<double> MinBnds;
    static std::vector<double> MaxBnds;

public:
    static void Run(std::vector<ExpData>* dataSet, int numberOfParameterGroups);
    static int GetNumberOfParameterGroups();

};



#endif // OPTIMIZATIONENGINE_HPP_INCLUDED

