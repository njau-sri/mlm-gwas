TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -llapack -lblas -lgfortran -lquadmath
QMAKE_CXXFLAGS += -std=c++11 -fopenmp
QMAKE_LFLAGS += -static -fopenmp

SOURCES += \
        main.cpp \
    emma.cpp \
    lapack.cpp \
    lsfit.cpp \
    statsutil.cpp \
    mlm_gwas.cpp \
    cmdline.cpp \
    pheno.cpp \
    vcf.cpp \
    util.cpp

HEADERS += \
    emma.h \
    lapack.h \
    lsfit.h \
    statsutil.h \
    cmdline.h \
    pheno.h \
    split.h \
    vcf.h \
    util.h
