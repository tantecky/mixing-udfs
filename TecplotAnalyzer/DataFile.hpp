#ifndef DATAFILE_HPP_INCLUDED
#define DATAFILE_HPP_INCLUDED

#include <map>
#include <string>

#include <fstream>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <limits>
#include <iomanip>

class DataFile
{
private:
    std::map<std::string, int> m_DataVariables;

    void LoadDigits(const char* const, std::ifstream&, std::string&);

public:
    void LoadDataFile(const char* const);

};

#endif // DATAFILE_HPP_INCLUDED
