#include <iostream>
#include <cstdlib>
// the general include for eo
#include <eo>
#include <es.h>

#include "ExpData.hpp"
#include "Integrator.hpp"
#include "ModEq.hpp"
#include "OptimizationEngine.hpp"


int main(int argc, char* argv[])
{
    try
    {
        double initalMass;
        int numberOfParameterGroups;

        if(argc != 4)
        {
            std::cerr << "Bad number of arguments" << std::endl;
            std::cerr << " expected: PyrolysisGA <file> <initialMass> <numberOfParameterGroups>" << std::endl;
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

        numberOfParameterGroups = std::atoi(argv[3]);

        if(numberOfParameterGroups < 1 || numberOfParameterGroups > 3)
        {
            std::cerr << "Wrong number: numberOfParameterGroups" << std::endl;

            return EXIT_FAILURE;
        }

        std::vector<ExpData> dataSet;

        ExpData::LoadExpData(argv[1], dataSet, initalMass); //initialMass

        OptimizationEngine::Run(&dataSet, numberOfParameterGroups);
    }
    catch(std::exception& e)
    {
        std::cerr << "Exception: " << std::endl;
        std::cerr << "    " << e.what() << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
