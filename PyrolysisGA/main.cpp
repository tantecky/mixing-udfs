#include <iostream>
#include <cstdlib>

#include "ExpData.hpp"


int main()
{
    try
    {
        std::vector<ExpData> dataSet;

        ExpData::LoadExpData("borovice_buk 5_15_25 vzduch.csv", dataSet);

        for(std::vector<ExpData>::const_iterator it = dataSet.begin(); it != dataSet.end(); ++it)
        {
            std::cout << *it;
        }
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
