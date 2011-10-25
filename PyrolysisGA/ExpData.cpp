#include "ExpData.hpp"

std::ostream& operator<<(std::ostream& stream, ExpData obj)
{
    return stream << "Time: " << obj.Time() << " Temp: " << obj.Temp() << " TG: " << obj.TG() << " HeatFlow: " << obj.HeatFlow() << std::endl;
}
