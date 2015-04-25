#ifndef MATRIX_H
#define MATRIX_H

#define ERR 1e-7

#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include "myexception.h"

template<typename T>
class Matrix
{
private:
    std::vector<T> mdata;
    std::vector<T*> data;

    int _row, _col;

    bool check(int i, int j) const {return i >=0 && i < _row && j >= 0 && j < _col;}
    bool checkr(int i) const {return i >=0 && i < _row;}
    bool checkc(int i) const {return i >=0 && i < _col;}

public:
    Matrix();
    Matrix(int row, int col = 1, const T& init_var = 0);
    Matrix(const Matrix<T>& mat);
    Matrix(const std::vector< Matrix<T> >& vecs);
    ~Matrix();

    //get / set
     int rows() const {return _row;}
     int cols() const {return _col;}
     int getDim() const {return _row;}
    void resize( int new_row,  int new_col);
    void swapRows(int ri, int rj);
    void swapCols(int ci, int cj);
    Matrix<T> getRow(int i);
    Matrix<T> getCol(int i);
    void setRow(int i, const T& var);
    void setCol(int i, const T& var);
    void setRow(int i, const Matrix<T>& vec);
    void setCol(int i, const Matrix<T>& vec);
    void fill(T var);
    std::string toString();

    //overload operator
    T& operator()(int row, int col = 0);
    const T& operator()(int row, int col = 0) const;

    Matrix<T>& operator=(const Matrix<T>& mat);

    Matrix<T> operator+(const Matrix<T>& mat) const;
    Matrix<T> operator+(const T& var) const;
    Matrix<T>& operator+=(const Matrix<T>& mat);
    Matrix<T>& operator+=(const T& var);

    Matrix<T>& operator-();
    Matrix<T> operator-(const Matrix<T>& mat) const;
    Matrix<T> operator-(const T& var) const;
    Matrix<T>& operator-=(const Matrix<T>& mat);
    Matrix<T>& operator-=(const T& var);

    Matrix<T> operator*(const Matrix<T>& mat) const;
    Matrix<T> operator*(const T& var) const;
    Matrix<T>& operator*=(const Matrix<T>& mat);
    Matrix<T>& operator*=(const T& var);

    Matrix<T> operator/(const T& var) const;
    Matrix<T>& operator/=(const T& var);

    bool operator==(const Matrix<T>& mat);
    bool operator!=(const Matrix<T>& mat) {return !(*this == mat);}
};

typedef Matrix<double> Matd;
typedef Matrix<double> Vecd;

#include "matrix.cpp"

#endif // MATRIX_H
