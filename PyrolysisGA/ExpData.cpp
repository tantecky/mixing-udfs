#include "ExpData.hpp"

double ExpData::m_AvgMassFrac = 0.0;
double ExpData::m_Beta = 5./60.;

const double ExpData::EqWeight[2][3][3]  =
{
    {
        //pine
        { 1.0},
        { 0.75, 0.25},      //cellulose + hemicelluloses, lignin
        { 0.43, 0.32, 0.25} //cellulose, hemicelluloses, lignin
    },
    {
        //oak
        { 1.0},
        { 0.75, 0.25},      //cellulose + hemicelluloses, lignin
        { 0.4, 0.35, 0.25}  //cellulose, hemicelluloses, lignin
    }

};

ExpData::Material ExpData::m_WoodType = ExpData::Pine;

double ExpData::AvgMassFrac()
{
    return ExpData::m_AvgMassFrac;
}

double ExpData::Beta()
{
    return ExpData::m_Beta;
}

const ExpData::Material& ExpData::WoodType()
{
    return m_WoodType;
}

void ExpData::LoadExpData(const char* const filename, std::vector<ExpData>& dataSet, double initialMass, Material woodType, double beta)
{
    m_WoodType = woodType;
    m_Beta = beta;

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
    double massFrac;

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

        massFrac = (initialMass-(initialMass*std::abs(data[2])/100.0))/initialMass;

        dataSet.push_back(ExpData(data[1], data[0] + 273.15, data[2], data[3], massFrac));

    }

    dataFile.close();

    double weightFrac = 0.0;
    double totTime = 0.0;

    for(unsigned int i = 0; i < dataSet.size(); i++)
    {
        totTime += dataSet[i].Time();
        weightFrac += dataSet[i].Time()*dataSet[i].MassFrac();
    }

    m_AvgMassFrac = weightFrac/totTime;

}

std::ostream& operator<<(std::ostream& stream, ExpData obj)
{
    return stream << "Time: " << obj.Time() << " Temp: " << obj.Temp() << " TG: " << obj.TG() << " HeatFlow: " << obj.HeatFlow() << " MassFrac: " << obj.MassFrac() <<std::endl;
}


