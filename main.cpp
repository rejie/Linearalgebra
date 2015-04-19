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
        Matrix<double> mat(50, 50);
        Vector<double> vec(50), res(50);

        for(int i=0; i<mat.rows(); ++i)
        {
            for(int j=0; j<mat.cols(); ++j)
            {
                cin >> c;
                mat(i, j) = c;
            }
        }
        cout << endl;
        for(int i=0; i<vec.getDim(); ++i)
        {
            cin >> c;
            vec(i) = c;
        }
        system("cls");
        res = solver<double>(mat, vec);
        cout << endl << res.toString() << endl;
    }
    catch(const MyException& e)
    {
        cout << e.what() << endl;
    }

    return 0;
}
