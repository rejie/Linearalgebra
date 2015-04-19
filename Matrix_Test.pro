#-------------------------------------------------
#
# Project created by QtCreator 2015-04-10T14:36:22
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = Matrix_Test
CONFIG   += console
CONFIG   -= app_bundle

DEFINES += _USE_MATH_DEFINES

TEMPLATE = app

INCLUDEPATH += $$PWD/Linearalgebra/

SOURCES += main.cpp \
    Linearalgebra/matrix.cpp \
    Linearalgebra/vector.cpp

HEADERS += Linearalgebra/matrix.h \
    Linearalgebra/myexception.h \
    solver.h \
    Linearalgebra/vector.h
