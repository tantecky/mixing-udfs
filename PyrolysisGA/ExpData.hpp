#ifndef EXPDATA_HPP_INCLUDED
#define EXPDATA_HPP_INCLUDED

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <cmath>

class ExpData
{
public:
    enum Material
    {
        Pine = 0,
        Oak = 1,
    };

private:
    static const int NumberOfRows = 4;

    double m_Time;
    double m_Temp;
    double m_TG;
    double m_HeatFlow;
    double m_MassFrac;
    static Material m_WoodType;

    static double m_AvgMassFrac;
    static double m_Beta;

public:

    const static double EqWeight[2][3][3];

    static double AvgMassFrac();
    static double Beta();
    static const Material& WoodType();

    static void LoadExpData(const char* const filename, std::vector<ExpData>& dataSet, double initialMass, Material woodType, double beta);

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

    inline const double& MassFrac()
    {
        return m_MassFrac;
    }

    ExpData(double time, double temp, double tg, double heatFlow, double massFrac)
        : m_Time(time), m_Temp(temp), m_TG(tg), m_HeatFlow(heatFlow), m_MassFrac(massFrac)
    {

    }

};

#endif // EXPDATA_HPP_INCLUDED
