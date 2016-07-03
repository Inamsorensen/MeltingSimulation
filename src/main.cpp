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

//  //Test conjugate gradient + linear system solve
//  std::vector<float> Amatrix;
//  std::vector<float> Bmatrix;
//  std::vector<float> xMatrix;
//  std::vector<float> x0Matrix;
//  float maxIter=4;
//  float minResidual=0.001;

//  Amatrix.push_back(2.0);
//  Amatrix.push_back(-1.0);
//  Amatrix.push_back(0.0);
//  Amatrix.push_back(-1.0);
//  Amatrix.push_back(2.0);
//  Amatrix.push_back(-1.0);
//  Amatrix.push_back(0.0);
//  Amatrix.push_back(-1.0);
//  Amatrix.push_back(2.0);

//  Bmatrix.push_back(1.0);
//  Bmatrix.push_back(1.0);
//  Bmatrix.push_back(1.0);

//  MathFunctions::conjugateGradient(&Amatrix, &Bmatrix, &x0Matrix, &xMatrix, maxIter, minResidual);

//  std::vector<float> xMatrix_linearSolver;

//  MathFunctions::linearSystemSolve(&Amatrix, &Bmatrix, &xMatrix_linearSolver);

//  std::cout<<"Solution "<<xMatrix[0]<<" "<<xMatrix[1]<<" "<<xMatrix[2]<<"\n";
//  std::cout<<"Solution "<<xMatrix_linearSolver[0]<<" "<<xMatrix_linearSolver[1]<<" "<<xMatrix_linearSolver[2]<<"\n";

  //Test singular value decomposition
//    std::vector<float> decomposeMatrix;
//    std::vector<float> U;
//    std::vector<float> S;
//    std::vector<float> V;

//    decomposeMatrix.push_back(2.031);
//    decomposeMatrix.push_back(0.162);
//    decomposeMatrix.push_back(2.021);
//    decomposeMatrix.push_back(0.162);
//    decomposeMatrix.push_back(1.245);
//    decomposeMatrix.push_back(0.0);
//    decomposeMatrix.push_back(2.021);
//    decomposeMatrix.push_back(0.0);
//    decomposeMatrix.push_back(2.561);

//    MathFunctions::singularValueDecomposition(&decomposeMatrix, &U, &S, &V);

//    std::cout<<"U: "<<U[0]<<" "<<U[1]<<" "<<U[2]<<"\n";
//    std::cout<<"   "<<U[3]<<" "<<U[4]<<" "<<U[5]<<"\n";
//    std::cout<<"   "<<U[6]<<" "<<U[7]<<" "<<U[8]<<"\n";

//    std::cout<<"S: "<<S[0]<<" "<<S[1]<<" "<<S[2]<<"\n";
//    std::cout<<"   "<<S[3]<<" "<<S[4]<<" "<<S[5]<<"\n";
//    std::cout<<"   "<<S[6]<<" "<<S[7]<<" "<<S[8]<<"\n";

//    std::cout<<"V: "<<V[0]<<" "<<V[1]<<" "<<V[2]<<"\n";
//    std::cout<<"   "<<V[3]<<" "<<V[4]<<" "<<V[5]<<"\n";
//    std::cout<<"   "<<V[6]<<" "<<V[7]<<" "<<V[8]<<"\n";

//  Eigen::Matrix3f decomposeMatrix;
//  decomposeMatrix(0,0)=1.0;
//  decomposeMatrix(0,1)=1.0;
//  decomposeMatrix(0,2)=2.0;
//  decomposeMatrix(1,0)=0.0;
//  decomposeMatrix(1,1)=0.0;
//  decomposeMatrix(1,2)=5.0;
//  decomposeMatrix(2,0)=6.0;
//  decomposeMatrix(2,1)=7.0;
//  decomposeMatrix(2,2)=8.0;

//  Eigen::Matrix3f U;
//  Eigen::Matrix3f X;
//  Eigen::Matrix3f V;

//  MathFunctions::singularValueDecomposition(&decomposeMatrix, &U, &X, &V);

////  std::cout<<"U:\n";
////  std::cout<<U<<"\n";
////  std::cout<<"X:\n";
////  std::cout<<X<<"\n";
////  std::cout<<"V:\n";
////  std::cout<<V<<"\n";

//  Eigen::Matrix3f R;
//  Eigen::Matrix3f S;

//  MathFunctions::polarDecomposition(&decomposeMatrix, &R, &S);

//  std::cout<<"R:\n";
//  std::cout<<R<<"\n";
//  std::cout<<"S:\n";
//  std::cout<<S<<"\n";

//  Eigen::Matrix3f V_conj;
//  Eigen::Matrix3f V_conjTrans;
//  V_conj=V.conjugate();
//  V_conjTrans=V_conj.transpose();

//  std::cout<<"V*:\n";
//  std::cout<<V_conjTrans<<"\n";


  //Test matrix constructors
//  Eigen::MatrixXf testMINRES_A(7,7);
//  testMINRES_A.setConstant(-1.0);

//  for (int i=0; i<7; i++)
//  {
//    testMINRES_A(i,i)=6.0;
//  }
//  Eigen::VectorXf testMINRES_b(7);
//  testMINRES_b.setOnes();

//  Eigen::MatrixXf testMINRES_A(3,3);
//  testMINRES_A(0,0)=1.0;
//  testMINRES_A(0,1)=2.0;
//  testMINRES_A(0,2)=4.0;
//  testMINRES_A(1,0)=2.0;
//  testMINRES_A(1,1)=3.0;
//  testMINRES_A(1,2)=2.0;
//  testMINRES_A(2,0)=4.0;
//  testMINRES_A(2,1)=2.0;
//  testMINRES_A(2,2)=6.0;

//  Eigen::VectorXf testMINRES_b(3);
//  testMINRES_b(0)=2.0;
//  testMINRES_b(1)=4.0;
//  testMINRES_b(2)=6.0;

//  Eigen::VectorXf x(3);
//  x.setZero();

//  Eigen::MatrixXf testMINRES_A(2,2);
//  testMINRES_A(0,0)=1.0;
//  testMINRES_A(0,1)=2.0;
//  testMINRES_A(1,0)=2.0;
//  testMINRES_A(1,1)=3.0;

//  Eigen::VectorXf testMINRES_b(2);
//  testMINRES_b(0)=1.0;
//  testMINRES_b(1)=2.0;

//  Eigen::VectorXf x(2);
//  x.setZero();

//  Eigen::MatrixXf* preconditioner=nullptr;

//  int maxLoops=10;

//  float tolerance=0.000000000001;

//  float shift=0.0;

//  bool show=true;



//  MathFunctions::MinRes(&testMINRES_A, &testMINRES_b, &x, preconditioner, shift, maxLoops, tolerance, show);

//  std::cout<<x<<"\n";

  return EXIT_SUCCESS;

}
