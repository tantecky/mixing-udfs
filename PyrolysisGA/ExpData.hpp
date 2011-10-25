#ifndef EXPDATA_HPP_INCLUDED
#define EXPDATA_HPP_INCLUDED

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <iostream>

class ExpData
{

private:
    static const int NumberOfRows = 4;

    double m_Time;
    double m_Temp;
    double m_TG;
    double m_HeatFlow;

public:
    static void LoadExpData(const char* const filename, std::vector<ExpData>& dataSet)
    {
        std::ifstream dataFile;

        dataFile.open(filename, std::ifstream::in);

        if(!dataFile.good())
        {
            std::string error("Unable to read file ");
            error += filename;

            throw std::invalid_argument(error);
        }

        std::string line;
        std::string cell;
        std::stringstream  lineStream;
        std::istringstream istring;

        double data[ExpData::NumberOfRows];

        while(getline(dataFile, line))
        {
            lineStream.clear();
            lineStream.str(line);

            for(int i = 0; i < ExpData::NumberOfRows; i++)
            {
                if(!getline(lineStream, cell, ','))
                    throw std::out_of_range("Bad number of rows");

                istring.clear();
                istring.str(cell);

                if(!(istring >> data[i]))
                {
                    std::string error("Unable convert ");
                    error += cell;
                    error += " to double";

                    throw std::invalid_argument(error);
                }

            }

            dataSet.push_back(ExpData(data[0], data[1], data[2], data[3]));

        }

        dataFile.close();
    }

    friend std::ostream& operator<<(std::ostream& stream, ExpData obj);

    inline const double& Time()
    {
        return m_Time;
    }

    inline const double& Temp()
    {
        return m_Temp;
    }

    inline const double& TG()
    {
        return m_TG;
    }

    inline const double& HeatFlow()
    {
        return m_HeatFlow;
    }

    ExpData(double time, double temp, double tg, double heatFlow)
        : m_Time(time), m_Temp(temp), m_TG(tg), m_HeatFlow(heatFlow)
    {

    }

};

#endif // EXPDATA_HPP_INCLUDED
