#include "OptimizationEngine.hpp"
#define STICK_BNDS

std::vector<ExpData>* OptimizationEngine::DataSet = (std::vector<ExpData>*)0;
double* OptimizationEngine::ModData = (double*)0;
double* OptimizationEngine::IntHeap = (double*)0;
int OptimizationEngine::NumberOfParameterGroups = 1;
std::vector<double> OptimizationEngine::MinBnds;
std::vector<double> OptimizationEngine::MaxBnds;


const double OptimizationEngine::PRE_FAC_MIN[3][3]  =
{
    { 1e14},             //NumberOfParameterGroups == 1
    { 1e5, 1e5},         //NumberOfParameterGroups == 2
    { 1e6, 1e6, 1e6}     //NumberOfParameterGroups == 3
};

const double OptimizationEngine::PRE_FAC_MAX[3][3]  =
{
    { 1e17},             //NumberOfParameterGroups == 1
    { 1e16, 1e16},       //NumberOfParameterGroups == 2
    { 1e20, 1e20, 1e20}  //NumberOfParameterGroups == 3
};

const double OptimizationEngine::ACT_ENER_MIN[3][3]  =
{
    { 1e5},              //NumberOfParameterGroups == 1
    { 1e5, 1e5},         //NumberOfParameterGroups == 2
    { 1e5, 1e5, 1e5}     //NumberOfParameterGroups == 3
};

const double OptimizationEngine::ACT_ENER_MAX[3][3]  =
{
    { 5e5},              //NumberOfParameterGroups == 1
    { 5e5, 5e5},         //NumberOfParameterGroups == 2
    { 5e5, 5e5, 5e5}     //NumberOfParameterGroups == 3
};

const double OptimizationEngine::MASS_FRAC_INF_MIN[3][3]  =
{
    { 0.001},            //NumberOfParameterGroups == 1
    { 0.001, 0.1},       //NumberOfParameterGroups == 2
    { 0.001, 0.001, 0.1} //NumberOfParameterGroups == 3
};

const double OptimizationEngine::MASS_FRAC_INF_MAX[3][3]  =
{
    { 0.35},             //NumberOfParameterGroups == 1
    { 0.05, 0.3},        //NumberOfParameterGroups == 2
    { 0.1, 0.1, 0.3}     //NumberOfParameterGroups == 3
};

const int OptimizationEngine::SEED = 1337;
const int OptimizationEngine::POP_SIZE = 250;                // Size of population
const int OptimizationEngine::MAX_GEN = 1200;                // Maximum number of generation before STOP
const unsigned int OptimizationEngine::PRINT_EVERY_SEC = 10; // print info to console every specified second
const double OptimizationEngine::HYPER_CUBE_RATE = 0.5;      // relative weight for hypercube Xover
const double OptimizationEngine::SEGMENT_RATE = 0.5;         // relative weight for segment Xover
const double OptimizationEngine::ALFA = 100.0;               // BLX coefficient
const double OptimizationEngine::EPSILON = 0.1;	             // range for real uniform mutation
const double OptimizationEngine::SIGMA = 0.3;	    	     // std dev. for normal mutation
const double OptimizationEngine::UNIFORM_MUT_RATE = 0.5;     // relative weight for uniform mutation
const double OptimizationEngine::DET_MUT_RATE = 0.5;         // relative weight for det-uniform mutation
const double OptimizationEngine::NORMAL_MUT_RATE = 0.5;      // relative weight for normal mutation
const double OptimizationEngine::P_CROSS = 0.8;	             // crossover probability
const double OptimizationEngine::P_MUT = 0.01;	             // mutation probability

int OptimizationEngine::GetNumberOfParameterGroups()
{
    return NumberOfParameterGroups;
}

void OptimizationEngine::Run(std::vector<ExpData>* dataSet, int numberOfParameterGroups)
{
    NumberOfParameterGroups = numberOfParameterGroups;
    DataSet = dataSet;
    ModData = new double[dataSet->size()];
    IntHeap = new double[dataSet->size()];

    std::vector<double> minBnds;

    for(int i = 0; i < NumberOfParameterGroups; i++)
    {
        minBnds.push_back(PRE_FAC_MIN[NumberOfParameterGroups-1][i]);        //A_min
        minBnds.push_back(ACT_ENER_MIN[NumberOfParameterGroups-1][i]);       //E_min
        minBnds.push_back(0.);                                               //NS_min
        minBnds.push_back(MASS_FRAC_INF_MIN[NumberOfParameterGroups-1][i]);  //yinf_min
    }


    std::vector<double> maxBnds;

    for(int i = 0; i < NumberOfParameterGroups; i++)
    {
        maxBnds.push_back(PRE_FAC_MAX[NumberOfParameterGroups-1][i]);     //A_max
        maxBnds.push_back(ACT_ENER_MAX[NumberOfParameterGroups-1][i]);     //E_max
        maxBnds.push_back(5.);                                     //NS_max
        maxBnds.push_back(MASS_FRAC_INF_MAX[NumberOfParameterGroups-1][i]);  //yinf_max
    }
    MinBnds = minBnds;
    MaxBnds = maxBnds;
    eoRealVectorBounds bnds(minBnds, maxBnds);
    //rng.reseed(SEED);
    eoEvalFuncPtr<Indi, double, const std::vector<double>& > plainEval(FitnessFce);
    eoEvalFuncCounter<Indi> eval(plainEval);

    eoPop<Indi> pop;

    InitPop(pop, plainEval);

    pop.sort();

    std::cout << "Initial Population:" << std::endl;
    std::cout << "----------------------------" << std::endl;
    std::cout << pop;
    //Selection operators

    // eoDetTournamentSelect<Indi> selectOne(30); //deterministic tournament selection
    // eoProportionalSelect<Indi> selectOne; //roulette wheel selection
    eoStochTournamentSelect<Indi>  selectOne(0.8); //Stochastic Tournament
    eoSelectPerc<Indi> select(selectOne, 2.0);
    /*It will select floor(rate*pop.size()) individuals and pushes them to
    the back of the destination population.*/
    // eoSelectPerc<Indi> select(selectOne); //rate == 1.0

    // Replacement
    //eoGenerationalReplacement<Indi> replace; //all offspring replace all parents
    eoPlusReplacement<Indi> replace; //the best from offspring+parents become the next generation
    // eoSSGAStochTournamentReplacement<Indi> replace(0.95);
    /*
    in which parents to be killed are chosen by a (reverse) stochastic tournament.
    Additional parameter (in the constructor) is the tournament rate, a double.
    */

    //Transformation
    eoSegmentCrossover<Indi> xoverS(bnds, ALFA);
    // uniform choice in hypercube built by the parents
    eoHypercubeCrossover<Indi> xoverA(bnds, ALFA);
    // Combine them with relative weights
    eoPropCombinedQuadOp<Indi> xover(xoverS, HYPER_CUBE_RATE);
    xover.add(xoverA, HYPER_CUBE_RATE);

    // Mutation
    // offspring(i) uniformly chosen in [parent(i)-epsilon, parent(i)+epsilon]
    eoUniformMutation<Indi>  mutationU(bnds, EPSILON);
    // k (=1) coordinates of parents are uniformly modified
    eoDetUniformMutation<Indi>  mutationD(bnds, EPSILON);
    // all coordinates of parents are normally modified (stDev SIGMA)
    double sigma = SIGMA;
    eoNormalMutation<Indi>  mutationN(bnds, sigma);
    // Combine them with relative weights
    eoPropCombinedMonOp<Indi> mutation(mutationU, UNIFORM_MUT_RATE);
    mutation.add(mutationD, DET_MUT_RATE);
    mutation.add(mutationN, NORMAL_MUT_RATE);

    eoSGATransform<Indi> transform(xover, P_CROSS, mutation, P_MUT);

    eoGenContinue<Indi> genCont(MAX_GEN);

    eoCombinedContinue<Indi> continuator(genCont);

    //statistics
    // but now you want to make many different things every generation
    // (e.g. statistics, plots, ...).
    // the class eoCheckPoint is dedicated to just that:
    // Declare a checkpoint (from a continuator: an eoCheckPoint
    // IS AN eoContinue and will be called in the loop of all algorithms)
    eoCheckPoint<Indi> checkpoint(continuator);

    // Create a counter parameter
    eoValueParam<unsigned> generationCounter(0, "Gen.");

    // Create an incrementor (sub-class of eoUpdater). Note that the
    // parameter's value is passed by reference,
    // so every time the incrementer is updated (every generation),
    // the data in generationCounter will change.
    eoIncrementor<unsigned> increment(generationCounter.value());
    // Add it to the checkpoint,
    // so the counter is updated (here, incremented) every generation
    checkpoint.add(increment);
    // now some statistics on the population:
    // Best fitness in population
    eoBestFitnessStat<Indi> bestStat;
    // Second moment stats: average and stdev
    eoSecondMomentStats<Indi> SecondStat;
    // Add them to the checkpoint to get them called at the appropriate time
    checkpoint.add(bestStat);
    checkpoint.add(SecondStat);
    // The Stdout monitor will print parameters to the screen ...
    eoStdoutMonitor monitor;

    // when called by the checkpoint (i.e. at every generation)
    // checkpoint.add(monitor);
    eoTimedMonitor timed(PRINT_EVERY_SEC);
    timed.add(monitor);
    checkpoint.add(timed);

    // the monitor will output a series of parameters: add them
    monitor.add(generationCounter);
    // monitor.add(eval); // because now eval is an eoEvalFuncCounter!
    monitor.add(bestStat);
    monitor.add(SecondStat);
    // A file monitor: will print parameters to ... a File, yes, you got it!
    eoFileMonitor fileMonitor("stats.xy", " ");

    // the checkpoint mechanism can handle multiple monitors
    checkpoint.add(fileMonitor);
    // the fileMonitor can monitor parameters, too, but you must tell it!
    fileMonitor.add(generationCounter);
    fileMonitor.add(bestStat);
    fileMonitor.add(SecondStat);

    //Final settings
    eoEasyEA<Indi> gga(checkpoint, eval, select, transform, replace);

    std::cout << "Working..." << std::endl;
    gga(pop); //GO!

    // Output
    // Print (sorted) intial population
    pop.sort();
    std::cout << "Final Population:" << std::endl;
    std::cout << "----------------------------" << std::endl;
    std::cout << pop;
    std::cout << "----------------------------" << std::endl;

    std::cout << "The Best member:" << std::endl;

    double fitness = pop[0].fitness();
    std::cout << "Fitness: " << fitness << std::endl;

    double* preFac = new double[NumberOfParameterGroups];
    double* actEner = new double[NumberOfParameterGroups];
    double* reacOrd = new double[NumberOfParameterGroups];
    double* massFracInf = new double[NumberOfParameterGroups];

    for(int i = 0; i < NumberOfParameterGroups; i++)
    {
        preFac[i] = (pop[0])[4*i];
        actEner[i] = (pop[0])[4*i+1];
        reacOrd[i] = (pop[0])[4*i+2];
        massFracInf[i] = (pop[0])[4*i+3];


        std::cout << "A" << i+1 << ": " << preFac[i] << std::endl;
        std::cout << "E" << i+1 << ": " << actEner[i] << std::endl;
        std::cout << "NS" << i+1 << ": " << reacOrd[i] << std::endl;
        std::cout << "yinf" << i+1 << ": " << massFracInf[i] << std::endl;
    }

    SaveResults(fitness, preFac, actEner, reacOrd, massFracInf);

    std::cout << "----------------------------" << std::endl;

    std::cout << "Results saved into results.xy" << std::endl;
    std::cout << "Statistics saved into stats.xy" << std::endl;

    delete[] ModData;
    delete[] IntHeap;

    delete[] preFac;
    delete[] actEner;
    delete[] reacOrd;
    delete[] massFracInf;
}

double OptimizationEngine::FitnessFce(const std::vector<double>& pars)
{
    double fitness;

    double deltaExpMod = 0.0;
    double deltaExpAvg = 0.0;
    double preFacSum = 0.0;

    ClearModData();

    for(int i = 0; i < NumberOfParameterGroups; i++)
    {
#ifdef STICK_BNDS
        for(int j = 0; j < 4; j++) //bounds checking
        {
            if(pars[4*i+j] > MaxBnds[4*i+j] ||  pars[4*i+j] < MinBnds[4*i+j])
                return 1e20;
        }
#endif

        double preFac = pars[4*i];
        double actEner = pars[4*i+1];
        double reacOrd = pars[4*i+2];
        double massFracInf = pars[4*i+3];



        preFacSum += (preFac - PRE_FAC_MIN[NumberOfParameterGroups-1][i])/PRE_FAC_MAX[NumberOfParameterGroups-1][i];

        ComputeModData(preFac, actEner, reacOrd, massFracInf, i);
    }

    for(unsigned int i = 1; i < DataSet->size(); i++)
    {
        deltaExpMod += pow((*DataSet)[i].MassFrac() - ModData[i], 2);
        deltaExpAvg += pow((*DataSet)[i].MassFrac() - ExpData::AvgMassFrac(), 2);
    }

    fitness = 0.985*(1.0 - (deltaExpAvg - deltaExpMod)/deltaExpAvg) + (0.015/NumberOfParameterGroups)*preFacSum;

    if(!std::isfinite(fitness) || fitness < 0.0)
        return 1e20;

    return fitness;
}

void OptimizationEngine::InitPop(eoPop<Indi>& pop, eoEvalFuncPtr<Indi, double, const std::vector<double>& > eval)
{
    for(int i = 0; i < POP_SIZE; i++)
    {
        Indi v;

        if(NumberOfParameterGroups == 1)
        {
            v.push_back(Logdist(rng, 17, false) + PRE_FAC_MIN[0][0]);                        //A
            v.push_back(Logdist(rng, 6, false) + 1e5);                                //E
            v.push_back(rng.uniform(0., 5.));                                         //NS
            v.push_back(rng.uniform(0., (*DataSet)[DataSet->size()-1].MassFrac() ));  //yinf
        }
        else if(NumberOfParameterGroups == 2)
        {
            v.push_back(Logdist(rng, 13, false) + PRE_FAC_MIN[1][0]);                        //A1
            v.push_back(Logdist(rng, 6, false) + 1e5);                                //E1
            v.push_back(rng.uniform(0., 5.));                                         //NS1
            v.push_back(rng.uniform(0., (*DataSet)[DataSet->size()-1].MassFrac() ));  //yinf1

            v.push_back(Logdist(rng, 7, false) + PRE_FAC_MIN[1][1]);                         //A2
            v.push_back(Logdist(rng, 6, false) + 1e5);                                //E2
            v.push_back(rng.uniform(0., 5.));                                         //NS2
            //v.push_back(rng.uniform(0., (*DataSet)[DataSet->size()-1].MassFrac())); //yinf2
            v.push_back(rng.uniform(0., 0.1));                                        //yinf2
        }
        else
        {
            v.push_back(Logdist(rng, 20, false) + PRE_FAC_MIN[1][0]);                        //A1
            v.push_back(Logdist(rng, 6, false) + 1e5);                                //E1
            v.push_back(rng.uniform(0., 5.));                                         //NS1
            v.push_back(rng.uniform(0., (*DataSet)[DataSet->size()-1].MassFrac() ));  //yinf1

            v.push_back(Logdist(rng, 15, false) + PRE_FAC_MIN[1][1]);                         //A2
            v.push_back(Logdist(rng, 6, false) + 1e5);                                 //E2
            v.push_back(rng.uniform(0., 5.));                                          //NS2
            v.push_back(rng.uniform(0., (*DataSet)[DataSet->size()-1].MassFrac()));    //yinf2

            v.push_back(Logdist(rng, 7, false) + PRE_FAC_MIN[1][1]);                          //A3
            v.push_back(Logdist(rng, 6, false) + 1e5);                                 //E3
            v.push_back(rng.uniform(0., 5.));                                          //NS3
            v.push_back(rng.uniform(0., (*DataSet)[DataSet->size()-1].MassFrac()));    //yinf3
        }

        eval(v);
        pop.push_back(v);
    }
}

void OptimizationEngine::SaveResults(double fitness, double* preFac, double* actEner, double* reacOrd, double* massFracInf)
{
    std::ofstream results;

    results.open("results.xy", std::ios::out | std::ios::trunc);

    if(!results.good())
    {
        throw new std::ios_base::failure("Unable to write results.xy");
    }

    results << "#Fitness:" << fitness << std::endl;

    for(int i = 0; i < NumberOfParameterGroups; i++)
    {
        results << "#A" << i+1 << ": " << preFac[i] << std::endl;
        results << "#E" << i+1 << ": " << actEner[i] << std::endl;
        results << "#NS" << i+1 << ": " << reacOrd[i] << std::endl;
        results << "#yinf" << i+1 << ": " << massFracInf[i] << std::endl;

        Integrator::Runge23(DataSet, IntHeap, preFac[i], actEner[i], reacOrd[i], massFracInf[i], i);

        std::ofstream line;

        std::stringstream fileName;

        fileName << "line" << (i+1) << ".xy";

        line.open(fileName.str().c_str(), std::ios::out | std::ios::trunc);

        line << "#Time Temp Exp Model"  << std::endl;

        for(unsigned int j = 0; j < DataSet->size(); j++)
        {
            line << (*DataSet)[j].Time() <<  " " << (*DataSet)[j].Temp() << " " << (*DataSet)[j].MassFrac() << " " << IntHeap[j] << std::endl;
        }

        line.close();

    }

    results << "#Time Temp Exp Model"  << std::endl;

    ClearModData();

    for(int i = 0; i < NumberOfParameterGroups; i++)
    {
        ComputeModData(preFac[i], actEner[i], reacOrd[i], massFracInf[i], i);
    }

    for(unsigned int i = 0; i < DataSet->size(); i++)
    {
        results << (*DataSet)[i].Time() <<  " " << (*DataSet)[i].Temp() << " " << (*DataSet)[i].MassFrac() << " " << ModData[i] << std::endl;
    }

    results.close();
}

double OptimizationEngine::Logdist(eoRng& gen, int maxExp, bool alsoNegative)
{

    double interval = 1. / (maxExp + 1.);

    double iniRand = gen.uniform();

    if(iniRand <= interval)
    {
        if(!alsoNegative)
            return gen.uniform();
        else
            return gen.flip(0.5) ? gen.uniform() : -gen.uniform();
    }

    for(int i = 1; i <= maxExp; i++)
    {
        if(iniRand <= interval*(i+1))
        {
            if(!alsoNegative)
                return gen.uniform(pow(10., i - 1), pow(10., i));
            else
                return gen.flip(0.5) ? gen.uniform(pow(10., i - 1), pow(10., i))
                       :
                       -gen.uniform(pow(10., i - 1), pow(10., i));
        }


    }

    if(!alsoNegative)
        return gen.uniform(pow(10., maxExp - 1), pow(10., maxExp));
    else
        return gen.flip(0.5) ? gen.uniform(pow(10., maxExp - 1), pow(10., maxExp))
               : -gen.uniform(pow(10., maxExp - 1), pow(10., maxExp));
}

inline void OptimizationEngine::ClearModData()
{
    std::memset(ModData, 0, sizeof(double) * DataSet->size());
    ModData[0] = (*DataSet)[0].MassFrac(); //initial condition
}

void OptimizationEngine::ComputeModData(double preFac, double actEner, double reacOrd, double massFracInf, int component)
{
    Integrator::Runge23(DataSet, IntHeap, preFac, actEner, reacOrd, massFracInf, component);

    for(unsigned int j = 1; j < DataSet->size(); j++)
    {
        ModData[j] += IntHeap[j];
    }
}

