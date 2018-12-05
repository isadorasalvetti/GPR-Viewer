
TEMPLATE = app
TARGET = Base

QT += gui opengl

CONFIG += c++11

INCLUDEPATH += .

# Input
HEADERS += glwidget.h mainwindow.h \
    trianglemesh.h \
    plyreader.h \
    vertexmatrix.h
FORMS += mainwindow.ui
SOURCES += glwidget.cpp main.cpp mainwindow.cpp \
    trianglemesh.cpp \
    plyreader.cpp \
    vertexmatrix.cpp

DISTFILES += \
    shaders/simpleshader.vert \
    shaders/simpleshader.frag

RESOURCES += \
    resources.qrc
