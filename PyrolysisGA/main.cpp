#include <iostream>
#include <cstdlib>


#include "ExpData.hpp"
#include "Integrator.hpp"
#include "ModEq.hpp"
#include "OptimizationEngine.hpp"


int main(int argc, char* argv[])
{
    try
    {
        double initalMass;

        if(argc != 3)
        {
            std::cerr << "Bad number of arguments" << std::endl;
            std::cerr << " expected: PyrolysisGA <file> <initialMass>" << std::endl;
            std::cerr << " file format: " << std::endl;
            std::cerr << "             columns: Temperature in Celsius,time,TG,heatFlow" << std::endl;

            return EXIT_FAILURE;
        }

        initalMass = std::atof(argv[2]);

        if(initalMass <= 0.0)
        {
            std::cerr << "Wrong number: initialMass" << std::endl;

            return EXIT_FAILURE;
        }

        std::vector<ExpData> dataSet;

        ExpData::LoadExpData(argv[1], dataSet, initalMass); //initialMass

        OptimizationEngine::Run(&dataSet);
    }
    catch(std::exception& e)
    {
        std::cerr << "Exception: " << std::endl;
        std::cerr << "    " << e.what() << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
