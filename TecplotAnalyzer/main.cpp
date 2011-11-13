#include <iostream>
#include <cstdlib>
#include <cstring>
#include "DataFile.hpp"
#include "DataPoint.hpp"
#include "DataSetInfo.hpp"

void WriteHelp()
{
    std::cout << "Invalid number of arguments" << std::endl;
    std::cout << "app inputFile outPutFile type=tecplot|post"<< std::endl;

    std::exit(EXIT_FAILURE);
}

int main(int argc, char* argv[])
{
    try
    {
        if(argc != 4)
        {
            WriteHelp();
        }

        DataFile dataFile;

        if(!std::strcmp(argv[3], "post"))
        {
            dataFile.LoadCFDPostDataFile(argv[1]);
            dataFile.WriteCFDPostInfoFile(argv[2]);
        }
        else if(!std::strcmp(argv[3], "tecplot"))
        {
            dataFile.LoadTecplotDataFile(argv[1]);
            dataFile.WriteTecplotInfoFile(argv[2]);
        }
        else
        {
            WriteHelp();
        }

        return EXIT_SUCCESS;

    }
    catch(std::exception& e)
    {

        std::cerr << "Exception has arised:" << std::endl;
        std::cerr << e.what() << std::endl;

        return EXIT_FAILURE;
    }


}

