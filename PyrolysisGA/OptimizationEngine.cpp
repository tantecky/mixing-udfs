#include "OptimizationEngine.hpp"

std::vector<ExpData>* OptimizationEngine::DataSet = (std::vector<ExpData>*)0;
double* OptimizationEngine::ModData = (double*)0;

const double OptimizationEngine::A_initial = 5.49e12;
const double OptimizationEngine::E_initial = 1.70e5;
const double OptimizationEngine::NS_initial = 3.56;
const double OptimizationEngine::yinf_initial = 0.231642;

const int OptimizationEngine::SEED = 1337;
const int OptimizationEngine::POP_SIZE = 100; // Size of population
const int OptimizationEngine::MAX_GEN = 5e4; // Maximum number of generation before STOP

const double OptimizationEngine::HYPER_CUBE_RATE = 0.5;     // relative weight for hypercube Xover
const double OptimizationEngine::SEGMENT_RATE = 0.5;  // relative weight for segment Xover
const double OptimizationEngine::ALFA = 10.0;     //BLX coefficient
const double OptimizationEngine::EPSILON = 0.1;	// range for real uniform mutation
const double OptimizationEngine::SIGMA = 0.3;	    	// std dev. for normal mutation
const double OptimizationEngine::UNIFORM_MUT_RATE = 0.5;  // relative weight for uniform mutation
const double OptimizationEngine::DET_MUT_RATE = 0.5;      // relative weight for det-uniform mutation
const double OptimizationEngine::NORMAL_MUT_RATE = 0.5;   // relative weight for normal mutation
const double OptimizationEngine::P_CROSS = 0.8;	// Crossover probability
const double OptimizationEngine::P_MUT = 0.5;	// mutation probability

void OptimizationEngine::Run(std::vector<ExpData>* dataSet)
{
    DataSet = dataSet;
    ModData = new double[dataSet->size()];

    rng.reseed(SEED);
    eoEvalFuncPtr<Indi, double, const std::vector<double>& > eval(FitnessFce);

    eoPop<Indi> pop;

    InitPop(pop, eval);

    pop.sort();

    std::cout << "Initial Population" << std::endl;
    std::cout << pop;

   // eoDetTournamentSelect<Indi> selectOne(30); //deterministic tournament selection
   // eoProportionalSelect<Indi> selectOne; //roulette wheel selection
    eoStochTournamentSelect<Indi>  selectOne(0.8); //Stochastic Tournament

    eoSelectPerc<Indi> select(selectOne, 2.0);
    /*It will select floor(rate*pop.size()) individuals and pushes them to
    the back of the destination population.*/
    //eoSelectPerc<Indi> select(selectOne); //rate == 1.0

    //eoGenerationalReplacement<Indi> replace; //all offspring replace all parents
    eoPlusReplacement<Indi> replace; //the best from offspring+parents become the next generation

    //eoSSGAStochTournamentReplacement<Indi> replace(0.8);
    /*
    in which parents to be killed are chosen by a (reverse) stochastic tournament.
    Additional parameter (in the constructor) is the tournament rate, a double.
    */

    //Transformation
    eoSegmentCrossover<Indi> xoverS(ALFA);
    // uniform choice in hypercube built by the parents
    eoHypercubeCrossover<Indi> xoverA(ALFA);
    // Combine them with relative weights
    eoPropCombinedQuadOp<Indi> xover(xoverS, HYPER_CUBE_RATE);
    xover.add(xoverA, HYPER_CUBE_RATE);

    // MUTATION
    // offspring(i) uniformly chosen in [parent(i)-epsilon, parent(i)+epsilon]
    eoUniformMutation<Indi>  mutationU(EPSILON);
    // k (=1) coordinates of parents are uniformly modified
    eoDetUniformMutation<Indi>  mutationD(EPSILON);
    // all coordinates of parents are normally modified (stDev SIGMA)
    double sigma = SIGMA;
    eoNormalMutation<Indi>  mutationN(sigma);
    // Combine them with relative weights
    eoPropCombinedMonOp<Indi> mutation(mutationU, UNIFORM_MUT_RATE);
    mutation.add(mutationD, DET_MUT_RATE);
    mutation.add(mutationN, NORMAL_MUT_RATE);

    eoSGATransform<Indi> transform(xover, P_CROSS, mutation, P_MUT);

    eoGenContinue<Indi> genCont(MAX_GEN);

    eoCombinedContinue<Indi> continuator(genCont);

    eoEasyEA<Indi> gga(continuator, eval, select, transform, replace);

    // Apply algo to pop - that's it!
    std::cout << "GO!" << std::endl;
    gga(pop);

    // OUTPUT
    // Print (sorted) intial population
    pop.sort();
    std::cout << "Final pop:" << std::endl;
    std::cout << pop << std::endl;

    std::cout << "The Best:" << std::endl;
    std::cout << "           A           E      NS      yinf" << std::endl;
    std::cout << pop[0] << std::endl;

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
        fitness -= delta*delta;
    }

    if(!std::isfinite(fitness))
        return -1e300;

    return fitness;
}

void OptimizationEngine::InitPop(eoPop<Indi>& pop, eoEvalFuncPtr<Indi, double, const std::vector<double>& > eval)
{
    for(int i = 0; i < POP_SIZE; i++)
    {
        Indi v;

        v.push_back(A_initial + rng.normal(A_initial*rng.uniform(0., 2.)));
        v.push_back(E_initial + rng.normal(E_initial*rng.uniform(0., 2.)));
        v.push_back(NS_initial + rng.normal(NS_initial*rng.uniform(0., 2.)));
        v.push_back(yinf_initial + rng.normal(yinf_initial*rng.uniform(0., 2.)));

        eval(v);
        pop.push_back(v);
    }
}



