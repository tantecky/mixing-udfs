#include "DataPoint.hpp"

DataPoint::DataPoint(double x, double y, double z,
                     double velLX, double velLY, double velLZ,
                     double velSX, double velSY, double velSZ,
                     double eps, double vos)
    : m_X(x), m_Y(y), m_Z(z),
    m_VelocityLiquidX(velLX), m_VelocityLiquidY(velLY), m_VelocityLiquidZ(velLZ),
    m_VelocitySolidX(velSX), m_VelocitySolidY(velSY), m_VelocitySolidZ(velSZ),
    m_Eps(eps), m_VOS(vos)
{

}

double DataPoint::SchillerNauman()
{
    double slip_x = m_VelocityLiquidX - m_VelocitySolidX;
    double slip_y = m_VelocityLiquidY - m_VelocitySolidY;
    double slip_z = m_VelocityLiquidZ - m_VelocitySolidZ;


    /*Euclidean norm of slip velocity*/
    double slip = sqrt(slip_x*slip_x + slip_y*slip_y + slip_z*slip_z);

    /*relative Reynolds number for the primary phase*/
    double reyp = RHO_L*slip*DIAMETER/MU_L;

    /*drag coefficient*/
    double cd0;

    /*Schiller-Nauman*/
    if(reyp > 1000.)
        cd0 = 0.44;
    else
        cd0=24.0*(1.0+0.15*pow(reyp, 0.687))/reyp;

    return cd0;
}

double DataPoint::Pinelli()
{
    /*---Pinelli Correction---*/

    /*energy dissipation rate*/
    double eps = m_Eps;

    /*Kolmogoroff microscale*/
    double kolscale = pow((NU_L*NU_L*NU_L)/eps, 0.25);

    double pinelli = 0.4*tanh(16.*kolscale/DIAMETER - 1.) + 0.6;

    double cd = SchillerNauman()/pow(pinelli, 2.0);

    return cd;
}

double DataPoint::Brucato()
{
    /*---Brucato Correction---*/

    /*energy dissipation rate*/
    double eps = m_Eps;

    /*Kolmogoroff microscale*/
    double kolscale = pow((NU_L*NU_L*NU_L)/eps, 0.25);

    double cd = SchillerNauman()*(1 + 8.76e-4*pow(DIAMETER/kolscale, 3.0));

    return cd;
}

double DataPoint::Khopkar()
{
    /*---Khopkar Correction---*/

    /*energy dissipation rate*/
    double eps = m_Eps;

    /*Kolmogoroff microscale*/
    double kolscale = pow((NU_L*NU_L*NU_L)/eps, 0.25);

    double cd = SchillerNauman()*(1 + 8.76e-5*pow(DIAMETER/kolscale, 3.0));

    return cd;
}

std::ostream& operator<<(std::ostream& stream, DataPoint obj)
{
   /* return stream <<
           obj.X() << "," << obj.Y() << "," << obj.Z() << "," <<
           obj.VelocityLiquidX() << "," << obj.VelocityLiquidY() << "," << obj.VelocityLiquidZ() << "," <<
           obj.VelocitySolidX() << "," << obj.VelocitySolidY() << "," << obj.VelocitySolidZ() << "," <<
           obj.Eps() << "," << obj.VOS()
           << std::endl;*/

    return stream <<
            obj.Y() << "," << obj.VOS() << "," <<
            obj.SchillerNauman() << "," << obj.Pinelli() << "," << obj.Brucato() << "," << obj.Khopkar()
            << std::endl;

}

