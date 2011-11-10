#include "DataFile.hpp"

void DataFile::LoadDataFile(const char* const filename)
{
    std::ifstream dataFile;

    dataFile.open(filename, std::ifstream::in);

    if(!dataFile.good())
    {
        std::string error("Unable to read file ");

        error += filename;
        error += " in function ";
        error += __FUNCTION__;
        error += " in line ";

        std::stringstream ss;
        ss << __LINE__;
        error += ss.str();

        throw std::invalid_argument(error);
    }

    std::string line;

    getline(dataFile, line);

    if(!dataFile.good() || line.find("VARIABLES =") == std::string::npos)
    {
        std::string error("Unable to find variables field ");

        error += filename;
        error += " in function ";
        error += __FUNCTION__;
        error += " in line ";

        std::stringstream ss;
        ss << __LINE__;
        error += ss.str();

        throw std::invalid_argument(error);
    }


    //VARIABLES = "
    int index = 0;
    line.replace(0, 13, "");
    line.erase(line.size() - 2, 2);

    //check unique?
    m_DataVariables.insert(std::pair<std::string, int>(line, index));

    // std::cout << std::hex <<  (int)line[0] << std::endl;

    // std::cout << std::hex <<  line << std::endl;


    //all other variables
    while(getline(dataFile, line) && line[0] == ' ')
    {
        line.erase(0, 2);
        line.erase(line.size() - 2, 2);

        index++;
        m_DataVariables.insert(std::pair<std::string, int>(line, index));

        //std::cout << std::hex <<  (int)line[1] << std::endl;

        // std::cout << std::hex <<  line << std::endl;

        // std::cout << std::hex <<  line.size() << std::endl;



        // std::cout << line << " " << line.size() <<  std::endl;
    }

    if(line.compare(0, 4, "ZONE") != 0)
    {
        std::string error("ZONE marker is missing ");

        error += filename;
        error += " in function ";
        error += __FUNCTION__;
        error += " in line ";

        std::stringstream ss;
        ss << __LINE__;
        error += ss.str();

        throw std::invalid_argument(error);
    }

    LoadDigits(filename, dataFile, line);
}

void DataFile::LoadDigits(const char* const filename, std::ifstream& dataFile, std::string& line)
{
    getline(dataFile, line);

    if(line.compare(0, 4, "DT=(") != 0)
    {
        std::string error("ZONE marker is missing ");

        error += filename;
        error += " in function ";
        error += __FUNCTION__;
        error += " in line ";

        std::stringstream ss;
        ss << __LINE__;
        error += ss.str();

        throw std::invalid_argument(error);
    }

    size_t pos = 0;
    int countOfColumns = 0;

    pos = line.find("DOUBLE", pos + 1);

    while(pos != std::string::npos)
    {
        countOfColumns++;
        pos = line.find("DOUBLE", pos + 1);
    }

    std::istringstream istring;
    size_t posSpaceStart;
    size_t posSpaceEnd;
    int countOfColumnsCheck;

    /*on the heap baby :) */
    double* rowOfDoubles = new double[countOfColumns];

    while(getline(dataFile, line))
    {

        posSpaceStart = 0;
        countOfColumnsCheck = 0;

        posSpaceStart = line.find_first_not_of(" \x0D", posSpaceStart);

        while(posSpaceStart != std::string::npos)
        {
            posSpaceEnd = line.find(' ', posSpaceStart + 1);

            //std::cout << line.substr(posSpaceStart, posSpaceEnd - posSpaceStart) << std::endl;

            std::istringstream istring;

            istring.clear();
            istring.str(line.substr(posSpaceStart, posSpaceEnd - posSpaceStart));

            if(!(istring >> rowOfDoubles[countOfColumnsCheck]))
            {
                std::string error("Unable convert to double");
                error += " in function ";
                error += __FUNCTION__;
                error += " in line ";

                std::stringstream ss;
                ss << __LINE__;
                error += ss.str();

                throw std::invalid_argument(error);
            }

            posSpaceStart = line.find_first_not_of(" \x0D", posSpaceEnd);
            countOfColumnsCheck++;

            if(countOfColumnsCheck  > countOfColumns)
            {
                std::string error("Error in columns arithmetics ");

                error += filename;
                error += " in function ";
                error += __FUNCTION__;
                error += " in line ";

                std::stringstream ss;
                ss << __LINE__;
                error += ss.str();

                delete[] rowOfDoubles;
                throw std::invalid_argument(error);
            }
        }

        //std::cout << rowOfDoubles[FindVariable(DataSetInfo::NameOfZ)] << std::endl;

        try
        {
            m_DataPoints.push_back(DataPoint(
                                       rowOfDoubles[FindVariable(DataSetInfo::NameOfX)], rowOfDoubles[FindVariable(DataSetInfo::NameOfY)], rowOfDoubles[FindVariable(DataSetInfo::NameOfZ)],
                                       rowOfDoubles[FindVariable(DataSetInfo::NameOfVelocityLiquidX)], rowOfDoubles[FindVariable(DataSetInfo::NameOfVelocityLiquidY)], rowOfDoubles[FindVariable(DataSetInfo::NameOfVelocityLiquidZ)],
                                       rowOfDoubles[FindVariable(DataSetInfo::NameOfVelocitySolidX)], rowOfDoubles[FindVariable(DataSetInfo::NameOfVelocitySolidY)], rowOfDoubles[FindVariable(DataSetInfo::NameOfVelocitySolidZ)],
                                       rowOfDoubles[FindVariable(DataSetInfo::NameOfEps)], rowOfDoubles[FindVariable(DataSetInfo::NameOfVOS)]
                                   ));
        }
        catch(std::invalid_argument& e)
        {
            delete[] rowOfDoubles;
            throw e;
        }

        //  break;
    }
    // std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << rowOfDoubles[8] << std::endl;

    delete[] rowOfDoubles;
    dataFile.close();

}

int DataFile::FindVariable(const std::string& var)
{
    m_Iter = m_DataVariables.find(var);

    if(m_Iter == m_DataVariables.end())
    {
        std::string error("Variable wasn't found ");

        error += var;
        error += " in function ";
        error += __FUNCTION__;
        error += " in line ";

        std::stringstream ss;
        ss << __LINE__;
        error += ss.str();

        throw std::invalid_argument(error);
    }

    return m_Iter->second;
}

void DataFile::WriteInfoFile(const char* const filename)
{
    std::ofstream outfile;
    outfile.open(filename);

    outfile << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    for(std::vector<DataPoint>::const_iterator it = m_DataPoints.begin(); it != m_DataPoints.end(); ++it)
    {
        outfile << *it;
    }

    outfile.close();
}
