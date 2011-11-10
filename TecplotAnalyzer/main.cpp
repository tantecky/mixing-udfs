#include <iostream>
#include <cstdlib>
#include "DataFile.hpp"
#include "DataPoint.hpp"
#include "DataSetInfo.hpp"

int main()
{
    try
    {
        DataFile dataFile;

        dataFile.LoadDataFile("Tomtom.dat");

        dataFile.WriteInfoFile("output.dat");

        DataFile dataFile2;

        dataFile2.LoadDataFile("Tomtom2.dat");

        dataFile2.WriteInfoFile("output2.dat");

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
