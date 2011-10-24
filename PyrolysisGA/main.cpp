#include <iostream>

using namespace std;


double suma(double a, double b)
{
    return a + b;
}

int main()
{
   // cout << "Hello world!" << endl;

    double neco = suma(5, 6);

    //neco = neco + 2;
    neco += 2;

    bool logic = false;

    if(logic == true)
    {
        cout << "ahoj";
    }
    else
    {
        cout << "cau";
    }



    return 0;
}
