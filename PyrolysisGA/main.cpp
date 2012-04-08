#include <iostream>
#include <cstdlib>
#include <eo>
#include <es.h>

#include "ExpData.hpp"
#include "Integrator.hpp"
#include "ModEq.hpp"

std::vector<ExpData> dataSet;

double A = 5.49e12;
double E = 1.70e5;
double NS = 3.56;
double yinf = 0.231642;

int main()
{
    try
    {
        ExpData::LoadExpData("borovice_buk 5,15,25 dusik.csv", dataSet, 5.0);

       /* for(std::vector<ExpData>::const_iterator it = dataSet.begin(); it != dataSet.end(); ++it)
        {
            std::cout << *it;
        }*/

        double* modData = new double[dataSet.size()];

        Integrator::Runge23(dataSet, modData, A, E, NS, yinf);

        for(unsigned int i = 0; i < dataSet.size(); i++)
        {
            std::cout << dataSet[i].Time() << "," << dataSet[i].MassFrac() << "," << modData[i] << std::endl;
        }

        delete[] modData;
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
