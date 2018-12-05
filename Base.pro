
TEMPLATE = app
TARGET = Base

QT += gui opengl

CONFIG += c++11

INCLUDEPATH += .
INCLUDEPATH += /home/isa/Documents/GP/GPR-Viewer/Eigen/

# Input
HEADERS += \
    glwidget.h \
    mainwindow.h \
    trianglemesh.h \
    plyreader.h \
    smoothingmatrix.h \
    ui_mainwindow.h

FORMS += mainwindow.ui

SOURCES += \
    glwidget.cpp \
    main.cpp \
    mainwindow.cpp \
    trianglemesh.cpp \
    plyreader.cpp \
    smoothingmatrix.cpp \

DISTFILES += \
    shaders/simpleshader.vert \
    shaders/simpleshader.frag \
    Eigen/CMakeLists.txt

RESOURCES += \
    resources.qrc
