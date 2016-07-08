#include <QtGui/QGuiApplication>
#include <iostream>

#include <math.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SVD>

#include "OpenGLWindow.h"
#include "ReadGeo.h"
#include "MathFunctions.h"

int main(int argc, char *argv[])
{
  QGuiApplication app(argc, argv);
  QSurfaceFormat format;
  format.setSamples(4);
  #if defined( DARWIN)
    format.setMajorVersion(3);
    format.setMinorVersion(2);
  #else
    format.setMajorVersion(4);
    format.setMinorVersion(3);
  #endif
  format.setProfile(QSurfaceFormat::CoreProfile);
  format.setDepthBufferSize(24);
  OpenGLWindow window;
  window.setFormat(format);
  std::cout<<"Profile is "<<format.majorVersion()<<" "<<format.minorVersion()<<"\n";
  window.resize(1024, 720);
  window.show();

  return app.exec();




  return EXIT_SUCCESS;

}
