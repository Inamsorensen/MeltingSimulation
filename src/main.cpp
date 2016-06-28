#include <QtGui/QGuiApplication>
#include <iostream>

#include <math.h>

#include "OpenGLWindow.h"
#include "ReadGeo.h"

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

//  ReadGeo* file=new ReadGeo("../HoudiniFiles/particles.geo");
//  std::vector<ngl::Vec3> positions;
//  std::vector<float> pointParameters;
//  file->getPointPositions(&positions);
//  file->getPointParameter_Float("temperature", &pointParameters);
////  file->getPointPositions(&positions);
//  float compLimit=file->getSimulationParameter_Float("gridSize");
//  std::cout<<compLimit<<"\n";
//  ngl::Vec3 gridOrigin=file->getSimulationParameter_Vec3("gridOrigin");
//  std::cout<<gridOrigin.m_x<<" "<<gridOrigin.m_y<<" "<<gridOrigin.m_z<<"\n";
//  delete file;

  return EXIT_SUCCESS;

}
