/*
 * Author: T. Antecky
 * Date: 28.10.2012
 * Description: Heat transfer coefficient sign changer for Fluent's profile file
                Obeying platform specific line endings and fixed/scientific float format
                I/O operation through stdin/stdout
 * Rev: 1.0
*/
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <stdint.h>

using namespace std;

static const char UnixLineEnd = '\n';
static char buffer[128];
static bool WinLineEnds = false;
static uint_fast64_t lineNum = 0;

static void DetectLineEnds()
{
    cin.getline(buffer, sizeof(buffer), UnixLineEnd);
    lineNum++;

    if(!cin.good() || cin.gcount() < 3)
    {
        cerr << "File is too short!" << endl;
        exit(EXIT_FAILURE);
    }

    cout << buffer << UnixLineEnd;

    if(buffer[cin.gcount()-2] == '\r')
        WinLineEnds = true;
}

static bool IsScientific()
{
    for(streamsize i = 0; i < cin.gcount(); i++)
    {
        if(buffer[i] == 'e' || buffer[i] == 'E')
            return true;
    }

    return false;
}

int main()
{
    char pattern1[] = "(heat-transfer-coef";
    char pattern2[] = ")";

    cout.precision(9);

    DetectLineEnds();

    while(true)
    {
        cin.getline(buffer, sizeof(buffer), UnixLineEnd);
        lineNum++;

        if(cin.good())
        {
            cout << buffer << UnixLineEnd;

            if(!strncmp(pattern1, buffer, sizeof(pattern1)-1))
            {

                while(true)
                {
                    cin.getline(buffer, sizeof(buffer), UnixLineEnd);
                    lineNum++;

                    if(cin.good())
                    {

                        if(!strncmp(pattern2, buffer, sizeof(pattern2)-1)) //closing bracket detected
                        {
                            cout << buffer << UnixLineEnd;
                            break;
                        }
                        else
                        {
                            char* pEnd;
                            double alfa = strtod(buffer,&pEnd);

                            if(buffer == pEnd || (*pEnd != '\0' /*LF*/ && *pEnd != '\r' /*CRLF*/) || errno == ERANGE || !isfinite(alfa))
                            {
                                cerr << "Fail to convert char into double (line: " << lineNum << ")" << endl;
                                return EXIT_FAILURE;
                            }
                            else
                            {
                                if(alfa != 0.0)
                                {
                                    if(IsScientific())
                                        cout << scientific << -1.0*alfa;
                                    else
                                        cout << fixed << -1.0*alfa;

                                    if(WinLineEnds)
                                        cout << '\r';

                                    cout << UnixLineEnd;
                                }

                                else
                                    cout << buffer << UnixLineEnd;
                            }
                        }

                    }
                    else //no closing bracket
                    {
                        cerr << "Closing tag \")\" was not found" << endl;
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

    cout.flush();

    return EXIT_SUCCESS;
}
