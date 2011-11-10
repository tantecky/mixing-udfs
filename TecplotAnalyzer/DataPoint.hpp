#ifndef DATAPOINT_HPP_INCLUDED
#define DATAPOINT_HPP_INCLUDED

#include <iostream>
#include <cmath>

#define DIAMETER 1.02e-3 /*diameter of solid particle*/
#define RHO_L 1011.44 /*density of liquid phase*/
#define MU_L 5e-3 /*dynamic viscosity of liquid phase*/
#define NU_L (MU_L/RHO_L) /*kinematic viscosity of liquid phase*/

class DataPoint
{
private:
    double m_X;
    double m_Y;
    double m_Z;

    double m_VelocityLiquidX;
    double m_VelocityLiquidY;
    double m_VelocityLiquidZ;

    double m_VelocitySolidX;
    double m_VelocitySolidY;
    double m_VelocitySolidZ;

    double m_Eps;
    double m_VOS;

public:

    inline const double& X()
    {
        return m_X;
    }

    inline const double& Y()
    {
        return m_Y;
    }

    inline const double& Z()
    {
        return m_Z;
    }

    inline const double& VelocityLiquidX()
    {
        return m_VelocityLiquidX;
    }

    inline const double& VelocityLiquidY()
    {
        return m_VelocityLiquidY;
    }

    inline const double& VelocityLiquidZ()
    {
        return m_VelocityLiquidZ;
    }


    inline const double& VelocitySolidX()
    {
        return m_VelocitySolidX;
    }

    inline const double& VelocitySolidY()
    {
        return m_VelocitySolidY;
    }

    inline const double& VelocitySolidZ()
    {
        return m_VelocitySolidZ;
    }

    inline const double& Eps()
    {
        return m_Eps;
    }

    inline const double& VOS()
    {
        return m_VOS;
    }

    double SchillerNauman();
    double Pinelli();
    double Brucato();
    double Khopkar();

    DataPoint(double x, double y, double z,
              double velLX, double velLY, double velLZ,
              double velSX, double velSY, double velSZ,
              double eps, double vos);

    friend std::ostream& operator<<(std::ostream& stream, DataPoint point);


};

#endif // DATAPOINT_HPP_INCLUDED


