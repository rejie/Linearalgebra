#include <iostream>
#include "myexception.h"
#include "vector.h"
#include "matrix.h"
#include "solver.h"
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[])
{

    try
    {
        double c;
        int size = 10;
        Matrix<double> mat(size, size);
        Matrix<double> v(size, size), d(size, size);


        for(int i=0; i<mat.rows(); ++i)
        {
            for(int j=0; j<mat.cols(); ++j)
            {
                cin >> c;
                mat(i, j) = c;
            }
        }
        system("cls");
        PM<double>(mat, v, d);
        cout << endl << v.toString() << endl;
        cout << endl << d.toString() << endl;

    }
    catch(const MyException& e)
    {
        cout << e.what() << endl;
    }

    return 0;
}
