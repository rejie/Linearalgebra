#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include "myexception.h"


template<typename T>
class Vector
{
private:
    std::vector<T> data;
    int dimensions ;
    bool check(int i) const {return i >=0 && i < dimensions;}
public:
    Vector();
    Vector(int d);
    Vector(int d, int init_var);
    ~Vector(){}

    //get data
    int getDim() const {return dimensions;}
    std::string toString();

    //overload operator
    T& operator()(int site);
    const T& operator() (int site) const;
    Vector<T>& operator=(const Vector<T>& vec);
    Vector<T> operator+(const Vector<T>& vec);
    Vector<T>& operator+=(const Vector<T>& vec);
    Vector<T>& operator-();
    Vector<T> operator-(const Vector<T>& vec);
    Vector<T>& operator-=(const Vector<T>& vec);
    T operator*(const Vector<T>& vec);
    Vector<T> operator*(const T& var);
    Vector<T>& operator*=(const T& var);
    Vector<T> operator/(const T& var);
    Vector<T>& operator/=(const T& var);
    bool operator==(const Vector<T>& vec);
    bool operator!=(const Vector<T>& vec);


};

#include "vector.cpp"

#endif // MATRIX_H
