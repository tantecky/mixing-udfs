#include <iostream>

// the general include for eo
#include <eo>
#include <ga.h>

typedef eoBit<double> Indi;

using namespace std;


const double VAHY[] = {1, 1, 2, 4, 12};
const double CENY[] = {1, 2, 2, 10, 4};


const int SEED = 42;	// seed for random number generator
const int T_SIZE = 3; // size for tournament selection
const int VEC_SIZE = 5; // Number of bits in genotypes
const int POP_SIZE = 10; // Size of population
const int MAX_GEN = 100; // Maximum number of generation before STOP

const double P_CROSS = 0.8;	// Crossover probability
const double P_MUT = 0.4;	// mutation probability

const double P_MUT_PER_BIT = 0.01;	// internal probability for bit-flip mutation

const double onePointRate = 0.5;     // rate for 1-pt Xover
const double twoPointsRate = 0.5;     // rate for 2-pt Xover
const double URate = 0.5;            // rate for Uniform Xover

const double bitFlipRate = 0.5;      // rate for bit-flip mutation
const double oneBitRate = 0.5;       // rate for one-bit mutation

double fitness(const vector<bool>& bagl)
{
    double vahaBatohu = bagl[0]*VAHY[0] + bagl[1]*VAHY[1] + bagl[2]*VAHY[2] + bagl[3]*VAHY[3] + bagl[4]*VAHY[4];
    double cenaBatohu = bagl[0]*CENY[0] + bagl[1]*CENY[1] + bagl[2]*CENY[2] + bagl[3]*CENY[3] + bagl[4]*CENY[4];

    if(vahaBatohu > 15.0)
        return -1e6;

    return cenaBatohu;
}

int mainn()
{
    try
    {

        rng.reseed(SEED);

        eoEvalFuncPtr<Indi, double, const vector<bool>& > eval(fitness);

        eoUniformGenerator<bool> uGen;
        eoInitFixedLength<Indi> random(VEC_SIZE, uGen);
        eoPop<Indi> pop(POP_SIZE, random);

        apply<Indi>(eval, pop);

        // OUTPUT
        // sort pop before printing it!
        pop.sort();
        // Print (sorted) intial population (raw printout)
        cout << "Initial Population" << endl;
        cout << pop;

        eoDetTournamentSelect<Indi> selectOne(T_SIZE);       // T_SIZE in [2,POP_SIZE]
        // is now encapsulated in a eoSelectPerc (entage)
        eoSelectPerc<Indi> select(selectOne, 3.0);// by default rate==1
        eoPlusReplacement<Indi> replace;

        // CROSSOVER
        // 1-point crossover for bitstring
        eo1PtBitXover<Indi> xover1;
        // uniform crossover for bitstring
        eoUBitXover<Indi> xoverU;
        // 2-pots xover
        eoNPtsBitXover<Indi> xover2(2);
        // Combine them with relative rates
        eoPropCombinedQuadOp<Indi> xover(xover1, onePointRate);
        xover.add(xoverU, URate);
        xover.add(xover2, twoPointsRate);

        // MUTATION
        // standard bit-flip mutation for bitstring
        eoBitMutation<Indi>  mutationBitFlip(P_MUT_PER_BIT);
        // mutate exactly 1 bit per individual
        eoDetBitFlip<Indi> mutationOneBit;
        // Combine them with relative rates
        eoPropCombinedMonOp<Indi> mutation(mutationBitFlip, bitFlipRate);
        mutation.add(mutationOneBit, oneBitRate);

        // The operators are  encapsulated into an eoTRansform object
        eoSGATransform<Indi> transform(xover, P_CROSS, mutation, P_MUT);

        eoGenContinue<Indi> genCont(MAX_GEN);
        eoCombinedContinue<Indi> continuator(genCont);

        eoEasyEA<Indi> gga(continuator, eval, select, transform, replace);

        cout << "\n        Here we go\n\n";
        gga(pop);

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

