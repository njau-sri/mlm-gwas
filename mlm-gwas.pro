TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += D:/boost_1_67_0

LIBS += -Ld:/mingw-w64-i686-lapack-3.8.0-3-any.pkg/mingw32/lib -llapack -lblas -lgfortran -lquadmath

QMAKE_LFLAGS += -static

SOURCES += \
        main.cpp \
    emma.cpp \
    lapack.cpp \
    lsfit.cpp \
    statsutil.cpp \
    mlm_gwas.cpp \
    cmdline.cpp \
    pheno.cpp \
    vcf.cpp

HEADERS += \
    emma.h \
    lapack.h \
    lsfit.h \
    statsutil.h \
    cmdline.h \
    pheno.h \
    split.h \
    vcf.h
