#ifndef DATASETINFO_HPP_INCLUDED
#define DATASETINFO_HPP_INCLUDED

#include <string>

class DataSetInfo
{
public:
    //Tecplot
    static const std::string NameOfX;
    static const std::string NameOfY;
    static const std::string NameOfZ;

    static const std::string NameOfVelocityLiquidX;
    static const std::string NameOfVelocityLiquidY;
    static const std::string NameOfVelocityLiquidZ;

    static const std::string NameOfVelocitySolidX;
    static const std::string NameOfVelocitySolidY;
    static const std::string NameOfVelocitySolidZ;

    static const std::string NameOfEps;
    static const std::string NameOfVOS;

    //CfdPost
    static const int NumberOfColumns;
    static const int NameOfXCFDPost;
    static const int NameOfYCFDPost;
    static const int NameOfZCFDPost;

    static const int NameOfVelocityLiquidXCFDPost;
    static const int NameOfVelocityLiquidYCFDPost;
    static const int NameOfVelocityLiquidZCFDPost;

    static const int NameOfVelocitySolidXCFDPost;
    static const int NameOfVelocitySolidYCFDPost;
    static const int NameOfVelocitySolidZCFDPost;

    static const int NameOfEpsCFDPost;
    static const int NameOfVOSCFDPost;


};



#endif // DATASETINFO_HPP_INCLUDED
