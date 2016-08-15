# This specifies the exe name
TARGET=MeltingSimulation
# where to put the .o files
OBJECTS_DIR=obj
# core Qt Libs to use add more here if needed.
QT+=gui opengl core

# as I want to support 4.8 and 5 this will set a flag for some of the mac stuff
# mainly in the types.h file for the setMacVisual which is native in Qt5
isEqual(QT_MAJOR_VERSION, 5) {
        cache()
        DEFINES +=QT5BUILD
}
# where to put moc auto generated files
MOC_DIR=moc
# on a mac we don't create a .app bundle file ( for ease of multiplatform use)
CONFIG-=app_bundle


# Auto include all .cpp files in the project src directory (can specifiy individually if required)
SOURCES+= $$PWD/src/main.cpp \
    src/SimulationController.cpp \
    src/Emitter.cpp \
    src/Particle.cpp \
    src/Grid.cpp \
    src/OpenGLWindow.cpp \
    src/MathFunctions.cpp \
    src/Grid_deviatoricVelocity.cpp \
    src/Grid_pressureVelocity.cpp \
    src/Grid_Temperature.cpp \
    src/ReadGeo.cpp \
    src/MinRes.cpp \
    src/Grid_updateParticleFromGrid.cpp \
    src/AlembicExport.cpp \
    src/Grid_interpolateParticleToGrid.cpp \
    src/Grid_deviatoricVelocity_New.cpp \
    src/Grid_interpolateGridToParticle.cpp

# same for the .h files
HEADERS+= $$PWD/include/ReadGeo.h \
    include/SimulationController.h \
    include/Emitter.h \
    include/Particle.h \
    include/Grid.h \
    include/OpenGLWindow.h \
    include/MathFunctions.h \
    include/CellCentre.h \
    include/CellFace.h \
    include/InterpolationData.h \
    include/AlembicExport.h


# and add the include dir into the search path for Qt and make
INCLUDEPATH +=./include
INCLUDEPATH +=/usr/local/include/eigen3/Eigen/
#INCLUDEPATH +=/usr/lib/gcc/x86_64-redhat-linux/4.8.2/include
INCLUDEPATH+=/usr/local/include/OpenEXR

linux*:INCLUDEPATH+= /public/devel/include
linux*:LIBS+=-L/public/devel/lib -lAlembic -lhdf5 -lhdf5_hl

linux*:LIBS+=-L/usr/local/lib -lHalf

macx:LIBS+=-lAlembic

LIBS+= -fopenmp

QMAKE_CXXFLAGS+= -fopenmp

# where our exe is going to live (root of project)
DESTDIR=./
# add the glsl shader files
OTHER_FILES+= shaders/*.glsl \
              shaders/*.vs \
              shaders/*.fs \
                                                        README.md
# were are going to default to a console app
CONFIG += console
# note each command you add needs a ; as it will be run as a single line
# first check if we are shadow building or not easiest way is to check out against current
!equals(PWD, $${OUT_PWD}){
        copydata.commands = echo "creating destination dirs" ;
        # now make a dir
        copydata.commands += mkdir -p $$OUT_PWD/shaders ;
        copydata.commands += echo "copying files" ;
        # then copy the files
        copydata.commands += $(COPY_DIR) $$PWD/shaders/* $$OUT_PWD/shaders/ ;
        # now make sure the first target is built before copy
        first.depends = $(first) copydata
        export(first.depends)
        export(copydata.commands)
        # now add it as an extra target
        QMAKE_EXTRA_TARGETS += first copydata
}
NGLPATH=$$(NGLDIR)
isEmpty(NGLPATH){ # note brace must be here
        message("including $HOME/NGL")
        include($(HOME)/NGL/UseNGL.pri)
}
else{ # note brace must be here
        message("Using custom NGL location")
        include($(NGLDIR)/UseNGL.pri)
}
