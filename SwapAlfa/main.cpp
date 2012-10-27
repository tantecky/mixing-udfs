#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <stdint.h>

using namespace std;

int main()
{
    char buffer[128];
    char pattern1[] = "(heat-transfer-coef";
    char pattern2[] = ")";
    uint_fast64_t lineNum = 0;

    cout.precision(9);

    while(true)
    {
        cin.getline(buffer, sizeof(buffer));
        lineNum++;

        if(cin.good())
        {
            cout << buffer << endl;

            if(!strncmp(pattern1, buffer, sizeof(pattern1)-1))
            {

                while(true)
                {
                    cin.getline(buffer, sizeof(buffer));
                    lineNum++;

                    if(cin.good())
                    {

                        if(!strncmp(pattern2, buffer, sizeof(pattern2)-1)) //closing bracket detected
                        {
                            cout << buffer << endl;
                            break;
                        }
                        else
                        {
                            char* pEnd;
                            double alfa = strtod(buffer,&pEnd);

                            if(buffer == pEnd || (*pEnd != '\0' /*LF*/ && *pEnd != '\r' /*CRLF*/) || errno == ERANGE || !isfinite(alfa))
                            {
                                cerr << "Fail convert char into double (line: " << lineNum << ")" << endl;
                                return EXIT_FAILURE;
                            }
                            else
                            {
                                if(alfa != 0.0)
                                    cout << -1.0*alfa << endl;
                                else
                                    cout << buffer << endl;
                            }
                        }

                    }
                    else //no closing bracket
                    {
                        cerr << "Closing tag \")\" was not found";
                        return EXIT_FAILURE;
                    }
                }

            }
        }
        else
        {
            break;
        }
    }

    return EXIT_SUCCESS;
}
