#include <iostream>
#include "myexception.h"
#include "matrix.h"
#include "solver.h"
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[])
{

    try
    {
       Vecd v(3);
       Matd
       double d;
       for(int i=0; i<v.getDim(); ++i)
       {
           cin >> d;
           v(i) = d;
       }

       cout << endl << (v * v).toString() << endl;
    }
    catch(const MyException& e)
    {
        cout << e.what() << endl;
    }

    return 0;
}
