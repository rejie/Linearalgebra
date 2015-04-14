#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <sstream>
#include <string>
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
    Matrix(int row, int col, const T& init_var);
    Matrix(int row, int col);
    Matrix(const Matrix<T>& mat);
    ~Matrix();

    //get data
     int rows() const {return _row;}
     int cols() const {return _col;}
    void resize( int new_row,  int new_col);
    void swapRows(int ri, int rj);
    void swapCols(int ci, int cj);
    std::string toString();

    //overload operator
    T& operator()(int row, int col);
    const T& operator()(int row, int col) const;

    Matrix<T>& operator=(const Matrix<T>& mat);

    Matrix<T> operator+(const Matrix<T>& mat);
    Matrix<T> operator+(const T& var);
    Matrix<T>& operator+=(const Matrix<T>& mat);
    Matrix<T>& operator+=(const T& var);

    Matrix<T>& operator-();
    Matrix<T> operator-(const Matrix<T>& mat);
    Matrix<T> operator-(const T& var);
    Matrix<T>& operator-=(const Matrix<T>& mat);
    Matrix<T>& operator-=(const T& var);

    Matrix<T> operator*(const Matrix<T>& mat);
    Matrix<T> operator*(const T& var);
    Matrix<T>& operator*=(const Matrix<T>& mat);
    Matrix<T>& operator*=(const T& var);

    Matrix<T> operator/(const T& var);
    Matrix<T>& operator/=(const T& var);

    bool operator==(const Matrix<T>& mat);
    bool operator!=(const Matrix<T>& mat) {return !(*this == mat);}
};

#include "matrix.cpp"

#endif // MATRIX_H
