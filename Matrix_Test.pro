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

TEMPLATE = app

INCLUDEPATH += $$PWD/Linearalgebra/

SOURCES += main.cpp

HEADERS += \
    matrix.h \
    myexception.h
