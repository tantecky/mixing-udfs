#include "OptimizationEngine.hpp"

std::vector<ExpData>* OptimizationEngine::DataSet = (std::vector<ExpData>*)0;
double* OptimizationEngine::ModData = (double*)0;

const double OptimizationEngine::A_initial = 5.49e12;
const double OptimizationEngine::E_initial = 1.70e5;
const double OptimizationEngine::NS_initial = 3.56;
const double OptimizationEngine::yinf_initial = 0.231642;

void OptimizationEngine::Run(std::vector<ExpData>* dataSet)
{
    DataSet = dataSet;
    ModData = new double[dataSet->size()];


    eoEvalFuncPtr<Indi, double, const std::vector<double>& > eval(FitnessFce);

    eoPop<Indi> pop;

    InitPop(pop, eval);

    pop.sort();

    std::cout << "Initial Population" << std::endl;
    std::cout << pop;

    delete[] ModData;

}

double OptimizationEngine::FitnessFce(const std::vector<double>& pars)
{
    double A = pars[0];
    double E = pars[1];
    double NS = pars[2];
    double yinf = pars[3];

    double fitness = 0.0;
    double delta;

    Integrator::Runge23(DataSet, ModData, A, E, NS, yinf);

    for(unsigned int i = 1; i < DataSet->size(); i++)
    {
        delta = (*DataSet)[i].MassFrac() - ModData[i];
        fitness += delta*delta;
    }

    return fitness;
}

void OptimizationEngine::InitPop(eoPop<Indi>& pop, eoEvalFuncPtr<Indi, double, const std::vector<double>& > eval)
{
    for(int i = 0; i < POP_SIZE; i++)
    {
        Indi v;

        //v.push_back(A_initial + rng.random(0, 1));
        double a = eo::rng.random<double>(0., 1.);
        /*v.push_back(E_initial + rng.random(0, 1)*E_initial*rng.normal());
        v.push_back(NS_initial + rng.random(0, 1)*NS_initial*rng.normal());
        v.push_back(yinf_initial + rng.random(0, 1)*yinf_initial*rng.normal());*/
        eval(v);
        pop.push_back(v);
    }
}
