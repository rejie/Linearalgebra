#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <sstream>
#include <string>
#include "myexception.h"


template<typename T>
class Vector
{
private:
    std::vector<T> data;
    int dimensions ;
    bool check(int i) const {return i >=0 && i < dimensions;}
public:
    Vector(int d);
    ~Vector(){};

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
    Vector<T>& operator-(const Vector<T>& vec);
    Vector<T>& operator-=(const Vector<T>& vec);
    T operator*(const Vector<T>& vec);
    Vector<T> operator*(const T& var);
    Vector<T&> operator*=(const T& var);
    Vector<T> operator/(const T& var);
    Vector<T>& operator/=(const T& var);
    bool operator==(const Vector<T>& vec);
    bool operator!=(const Vector<T>& vec);

    //vector operation
    T norm(const Vector<T>& vec);
    Vector<T> normal(const Vector<T>& vec);
    Vector<T> cross(const Vector<T>& vec1 , const Vector<T>& vec2);
    T com(const Vector<T>& vec1 , const Vector<T>& vec2);
    Vector<T> (const Vector<T>& vec1 , const Vector<T>& vec2);
    T area(const Vector<T>& vec1 , const Vector<T>& vec2);
    bool isParallel(const Vector<T>& vec1 , const Vector<T>& vec2);
    bool isOrthogonal(const Vector<T>& vec1 , const Vector<T>& vec2);
    double angle(const Vector<T>& vec1 , const Vector<T>& vec2);
    Vector<T> pn(const Vector<T>& vec1 , const Vector<T>& vec2);
    std::vector< Vector<T> > ob(const std::vector< Vector<T> > v);
    Vector<T> proj(const Vector<T>& vec1 , const Vector<T>& vec2);

};

#include "vector.cpp"

#endif // MATRIX_H
