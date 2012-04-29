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
        ExpData::Material woodType;
        double beta;

        if(argc != 6)
        {
            std::cerr << "Bad number of arguments" << std::endl;
            std::cerr << " expected: PyrolysisGA <file> <initialMass> <numberOfParameterGroups> <woodType=pine|oak> <beta>" << std::endl;
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

        if(!std::strcmp(argv[4], "pine"))
        {
            woodType = ExpData::Pine;
        }
        else if(!std::strcmp(argv[4], "oak"))
        {
            woodType = ExpData::Oak;
        }
        else
        {
            std::cerr << "Unkown wood type" << std::endl;

            return EXIT_FAILURE;
        }

        beta = std::atof(argv[5]);

        if(initalMass <= 0.0)
        {
            std::cerr << "Wrong number: beta" << std::endl;

            return EXIT_FAILURE;
        }

        std::vector<ExpData> dataSet;

        ExpData::LoadExpData(argv[1], dataSet, initalMass, woodType, beta); //initialMass

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
