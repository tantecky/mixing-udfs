#ifndef DATAFILE_HPP_INCLUDED
#define DATAFILE_HPP_INCLUDED

#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <limits>
#include <iomanip>

#include "DataSetInfo.hpp"
#include "DataPoint.hpp"


class DataFile
{
private:
    std::map<std::string, int> m_DataVariables;
    std::map<std::string, int>::const_iterator m_Iter;
    std::vector<DataPoint> m_DataPoints;

    void LoadTecplotDigits(const char* const, std::ifstream&, std::string&);
    int FindVariable(const std::string&);


public:
    void LoadTecplotDataFile(const char* const);
    void LoadCFDPostDataFile(const char* const);
    void WriteTecplotInfoFile(const char* const);
    void WriteCFDPostInfoFile(const char* const);

};

#endif // DATAFILE_HPP_INCLUDED
