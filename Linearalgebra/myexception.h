#ifndef MYEXCEPTION_H
#define MYEXCEPTION_H

#include <stdexcept>

class MyException : public std::runtime_error
{
public:
    MyException(const char* msg) : std::runtime_error(msg) {}
};

#endif // EXCEPTION_H
