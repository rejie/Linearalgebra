#include <iostream>
#include "myexception.h"
#include "matrix.h"

using namespace std;

int main(int argc, char *argv[])
{

    try
    {
        Matrix<int> mat = identity<int>(10);

        mat.swapRows(0, 2);
        mat.swapCols(0, 1);
        mat = transpos(mat);

        cout << mat.toString();
    }
    catch(const MyException& e)
    {
        cout << e.what() << endl;
    }

    return 0;
}
