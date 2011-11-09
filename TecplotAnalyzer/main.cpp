#include <iostream>
#include <cstdlib>
#include "DataFile.hpp"

int main()
{
    try
    {
        DataFile dataFile;

        dataFile.LoadDataFile("Tomtom2.dat");


        return 0;
    }
    catch(std::exception& e)
    {

        std::cerr << "Exception has arised:" << std::endl;
        std::cerr << e.what() << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
