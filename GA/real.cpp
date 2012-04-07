#include <iostream>

// the general include for eo
#include <eo>
#include <es.h>

using namespace std;

typedef eoReal<eoMinimizingFitness> Indi;

const int SEED = 42;	// seed for random number generator
const int T_SIZE = 10; // size for tournament selection
const int VEC_SIZE = 2; // Number of bits in genotypes
const int POP_SIZE = 50; // Size of population
const int MAX_GEN = 50; // Maximum number of generation before STOP

// some parameters for chosing among different operators
const double hypercubeRate = 0.5;     // relative weight for hypercube Xover
const double segmentRate = 0.5;  // relative weight for segment Xover
const double uniformMutRate = 0.5;  // relative weight for uniform mutation
const double detMutRate = 0.5;      // relative weight for det-uniform mutation
const double normalMutRate = 0.5;   // relative weight for normal mutation
const double EPSILON = 0.01;	// range for real uniform mutation
double SIGMA = 0.3;	    	// std dev. for normal mutation

const float P_CROSS = 0.8;	// Crossover probability
const float P_MUT = 0.1;	// mutation probability

/*const double n = 0.22513691044544;
const double m = 42125.47;*/

const double SHEAR[] =
{
    48.6,
    73,
    97.3,
    145.9,
    206.7,
    304,
    437.8,
    632.4,
    899.9,
    1301.3
};

const double VISC[] =
{
    2072.9,
    1522.2,
    1215.1,
    890.8,
    674.3,
    499.6,
    378.5,
    282.2,
    211.8,
    162.5
};

double logdist(eoRng& gen, int maxExp, bool alsoNegative)
{
    /* if(maxExp < 1)
         throw new */

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

    //unlikely
    if(!alsoNegative)
        return gen.uniform(pow(10., maxExp - 1), pow(10., maxExp));
    else
        return gen.flip(0.5) ? gen.uniform(pow(10., maxExp - 1), pow(10., maxExp))
               : -gen.uniform(pow(10., maxExp - 1), pow(10., maxExp));
}

double fitness(const vector<double>& pars)
{
    double sum = 0;

    for(int i = 0; i < 10; i++)
    {
        sum += pow(  pars[0]*pow(SHEAR[i], pars[1] - 1.0) - VISC[i], 2.0);

    }

    return sum;

}

void initPop(eoPop<Indi>& pop, eoEvalFuncPtr<Indi, double, const vector<double>& >& eval)
{

}

int main()
{
    try
    {
        rng.reseed(SEED);
        eoEvalFuncPtr<Indi, double, const vector<double>& > eval(  fitness );


        eoPop<Indi> pop;
        // initPop(pop, eval);


        for(int i = 0; i < POP_SIZE; i++)
        {
            Indi v;

            v.push_back(logdist(rng, 8, false));
            v.push_back(rng.uniform(0, 1));

            eval(v);
            pop.push_back(v);
        }

        // OUTPUT
        // sort pop before printing it!
        pop.sort();
        // Print (sorted) intial population (raw printout)
        cout << "Initial Population" << endl;
        cout << pop;

        eoDetTournamentSelect<Indi> selectOne(T_SIZE);

        eoSelectPerc<Indi> select(selectOne, 3.0);// by default rate==1
        eoPlusReplacement<Indi> replace;

        vector<double> min;
        min.push_back(0.);
        min.push_back(0.);

        vector<double> max;
        max.push_back(1e8);
        max.push_back(1.);


        eoRealVectorBounds bnds(min, max);

        eoSegmentCrossover<Indi> xoverS(bnds, 10);
        // uniform choice in hypercube built by the parents
        eoHypercubeCrossover<Indi> xoverA(bnds, 10);
        // Combine them with relative weights
        eoPropCombinedQuadOp<Indi> xover(xoverS, segmentRate);
        xover.add(xoverA, hypercubeRate);


// MUTATION
        // offspring(i) uniformly chosen in [parent(i)-epsilon, parent(i)+epsilon]
        eoUniformMutation<Indi>  mutationU(bnds, EPSILON);
        // k (=1) coordinates of parents are uniformly modified
        eoDetUniformMutation<Indi>  mutationD(bnds, EPSILON);
        // all coordinates of parents are normally modified (stDev SIGMA)
        eoNormalMutation<Indi>  mutationN(bnds, SIGMA);
        // Combine them with relative weights
        eoPropCombinedMonOp<Indi> mutation(mutationU, uniformMutRate);
        mutation.add(mutationD, detMutRate);
        mutation.add(mutationN, normalMutRate);


        eoSGATransform<Indi> transform(xover, P_CROSS, mutation, P_MUT);

        eoGenContinue<Indi> genCont(MAX_GEN);

        eoCombinedContinue<Indi> continuator(genCont);

        eoEasyEA<Indi> gga(continuator, eval, select, transform, replace);

        // Apply algo to pop - that's it!
        cout << "\n        Here we go\n\n";
        gga(pop);

// OUTPUT
        // Print (sorted) intial population
        pop.sort();
        cout << "FINAL Population\n" << pop << endl;

        return 0;
    }
    catch(exception& e)
    {
        cerr << "Exception: " << e.what() << endl;
        return 1;
    }
}




