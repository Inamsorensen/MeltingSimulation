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

//  Eigen::SparseMatrix<double> A(2,2);
//  Eigen::VectorXd b(2);
//  Eigen::MatrixXf preconditioner;
//  Eigen::VectorXd x(2);
//  x.setZero();

//  A.insert(0,0)=1.0;
//  A.insert(1,0)=2.0;
//  A.insert(0,1)=2.0;
//  A.insert(1,1)=3.0;

//  b(0)=1.0;
//  b(1)=2.0;

//  float shift=0.0;
//  int maxLoops=100;
//  float tolerance=0.000000001;
//  bool show=true;
////  MathFunctions::MinRes(A, b, x, preconditioner, shift, maxLoops, tolerance, show);
//  MathFunctions::conjugateGradient(A, b, x, maxLoops, tolerance);
//  std::cout<<A<<"\n";
//  std::cout<<b<<"\n";
//  std::cout<<x<<"\n";

//  Eigen::Matrix3f A;
//  Eigen::Matrix3f B;

//  A=A.Identity();

//  A(0,0)=1.0;
//  A(1,0)=2.0;
//  A(2,0)=4.0;
//  A(0,1)=2.0;
//  A(0,2)=4.0;
//  A(1,1)=3.0;
//  A(1,2)=2.0;
//  A(2,1)=2.0;
//  A(2,2)=6.0;

//  B(0,0)=1.0;
//  B(1,0)=2.0;
//  B(2,0)=3.0;
//  B(0,1)=1.0;
//  B(0,2)=2.0;
//  B(1,1)=3.0;
//  B(1,2)=1.0;
//  B(2,1)=2.0;
//  B(2,2)=3.0;

//  Eigen::Matrix3f C=MathFunctions::matrixElementMultiplication(A, B);

//  Eigen::Matrix3f C_trans=C.transpose();
//  Eigen::Matrix3f C_transInv=C_trans.inverse();
//  Eigen::Matrix3f C_inv=C.inverse();
//  Eigen::Matrix3f C_invTrans=C_inv.transpose();

//  std::cout<<A<<"\n";
//  std::cout<<"----------\n";
//  std::cout<<B<<"\n";
//  std::cout<<"----------\n";
//  std::cout<<C<<"\n";
//  std::cout<<"----------\n";
//  std::cout<<C_transInv<<"\n";
//  std::cout<<"----------\n";
//  std::cout<<C_invTrans<<"\n";
//  std::cout<<"----------\n";



  return EXIT_SUCCESS;

}
