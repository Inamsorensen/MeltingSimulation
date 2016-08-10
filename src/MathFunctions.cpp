#include "MathFunctions.h"

#include <stdexcept>
#include <cmath>
#include <math.h>
#include <iostream>
#include <limits>

//----------------------------------------------------------------------------------------------------------------------

int MathFunctions::getVectorIndex(int i, int j, int k, int _noCells)
{
  int vectorIndex=i+(_noCells*j)+(_noCells*_noCells*k);
  return vectorIndex;
}

//----------------------------------------------------------------------------------------------------------------------

Eigen::Vector3i MathFunctions::getParticleGridCell(Eigen::Vector3f _particlePosition, float _cellSize, Eigen::Vector3f _gridEdgeOrigin)
{
  Eigen::Vector3i index;

  //Find grid indices from particle position.
  Eigen::Vector3f indexParticle=(1.0/_cellSize)*(_particlePosition-_gridEdgeOrigin);

  //Find which cell the particle is in
  index(0)=floor(indexParticle(0));
  index(1)=floor(indexParticle(1));
  index(2)=floor(indexParticle(2));

  return index;
}

//----------------------------------------------------------------------------------------------------------------------

float MathFunctions::calcCubicBSpline(float _x)
{
  /*Cubic B-Spline
  --------------------------------------------------------------------------------------------------
  N(x) = (|x|^3)/2 - x^2 + 2/3            if 0<=|x|<1
       = -(|x|^3)/6 + x^2 - 2*|x| + 4/3   if 1<=|x|<2
       = 0                              otherwise
  --------------------------------------------------------------------------------------------------
  */

  float result=0.0;
  float absX=std::abs(_x);

  if (absX<1.0)
  {
    result=0.5*(pow(absX,3));
    result-=pow(_x,2);
    result+=(2.0/3.0);
  }

  else if (absX>=1.0 && absX<2.0)
  {
    result=(-1.0/6.0)*(pow(absX,3));
    result+=pow(_x,2);
    result-=2.0*absX;
    result+=(4.0/3.0);
  }

  return result;
}

//----------------------------------------------------------------------------------------------------------------------

float MathFunctions::calcCubicBSpline_Diff(float _x)
{
  /*Cubic B-Spline differentiated
  --------------------------------------------------------------------------------------------------
  N(x) = (3/2)*(|x|^2)*signFunction(x) - 2x                       if 0<=|x|<1
       = -(1/2)*(|x|^2)*signFunction(x) + 2x - 2*signFunction(x)  if 1<=|x|<2
       = 0                                                        otherwise
  --------------------------------------------------------------------------------------------------
  */

  float result=0.0;
  float absX=std::abs(_x);

  float signX=MathFunctions::signFunction(_x);

  if (absX<1.0)
  {
    result=(3.0/2.0)*(pow(absX,2))*signX;
    result-=(2.0*_x);
  }

  else if (absX>=1.0 && absX<2.0)
  {
    result=(-0.5)*(pow(absX,2))*signX;
    result+=(2.0*_x);
    result-=(2.0*signX);
  }

  return result;
}

//----------------------------------------------------------------------------------------------------------------------

float MathFunctions::calcCubicBSpline_Integ(int _faceDirection, int _iIndexIncrement, int _jIndexIncrement, int _kIndexIncrement)
{
  /*Cubic B Spline Integrated
  --------------------------------------------------------------------------------------------------
  integrate(N(x)) = (1/8)*(|x|^4)*sign(x) - (1/3)*(x^3) + (2/3)*x                       if 0<=|x|<1
                  = -(1/24)*(|x|^4)*sign(x) + (1/3)*(x^3) - (|x|^2)*sign(x) + (4/3)*x   if 1<=|x|<2
                  = 0                                                                   otherwise

  Gives tabulated results as follows
  x1=1    x0=0    :  11/24         :  0.458333
  x1=0.5  x0=0    :  115/384       :  0.299479
  x1=1    x0=0.5  :  61/384        :  0.158854
  x1=2    x0=1    :  1/24          :  0.041667
  x1=1.5  x0=1    :  5/128         :  0.039063
  x1=2    x0=1.5  :  2.6048*10^-3  :  0.002605

  Which result is used depends on the i,j,k increments and which face is being calculated for
  --------------------------------------------------------------------------------------------------
  */

  float result=0.0;

  float fullIntegralNear=0.458333;
  float halfLowerIntegralNear=0.299479;
  float halfUpperIntegralNear=0.158854;
  float fullIntegralFar=0.041667;
  float halfLowerIntegralFar=0.039063;
  float halfUpperIntegralFar=0.002605;

  float twiceHalfLowerNear=2.0*halfLowerIntegralNear;
  float acrossNearFar=halfUpperIntegralNear + halfLowerIntegralFar;


  float integralX=0.0;
  float integralY=0.0;
  float integralZ=0.0;

  //Determine face
  switch (_faceDirection) {
  case 0:
  {
    switch (_iIndexIncrement)
    {
    case 0:
      integralX=fullIntegralNear;
      break;
    case 1:
      integralX=fullIntegralFar;
      break;
    case 2:
      integralX=0.0;
      break;
    case -1:
      integralX=fullIntegralNear;
      break;
    case -2:
      integralX=fullIntegralFar;
      break;
    default:
      break;
    }

    switch (_jIndexIncrement)
    {
    case 0:
      integralY=twiceHalfLowerNear;
      break;
    case 1:
      integralY=acrossNearFar;
      break;
    case 2:
      integralY=halfUpperIntegralFar;
      break;
    case -1:
      integralY=acrossNearFar;
      break;
    case -2:
      integralY=halfUpperIntegralFar;
      break;
    default:
      break;
    }

    switch (_kIndexIncrement)
    {
    case 0:
      integralZ=twiceHalfLowerNear;
      break;
    case 1:
      integralZ=acrossNearFar;
      break;
    case 2:
      integralZ=halfUpperIntegralFar;
      break;
    case -1:
      integralZ=acrossNearFar;
      break;
    case -2:
      integralZ=halfUpperIntegralFar;
      break;
    default:
      break;
    }

    break;
  }
  case 1:
  {
    switch (_jIndexIncrement)
    {
    case 0:
      integralY=fullIntegralNear;
      break;
    case 1:
      integralY=fullIntegralFar;
      break;
    case 2:
      integralY=0.0;
      break;
    case -1:
      integralY=fullIntegralNear;
      break;
    case -2:
      integralY=fullIntegralFar;
      break;
    default:
      break;
    }

    switch (_iIndexIncrement)
    {
    case 0:
      integralX=twiceHalfLowerNear;
      break;
    case 1:
      integralX=acrossNearFar;
      break;
    case 2:
      integralX=halfUpperIntegralFar;
      break;
    case -1:
      integralX=acrossNearFar;
      break;
    case -2:
      integralX=halfUpperIntegralFar;
      break;
    default:
      break;
    }

    switch (_kIndexIncrement)
    {
    case 0:
      integralZ=twiceHalfLowerNear;
      break;
    case 1:
      integralZ=acrossNearFar;
      break;
    case 2:
      integralZ=halfUpperIntegralFar;
      break;
    case -1:
      integralZ=acrossNearFar;
      break;
    case -2:
      integralZ=halfUpperIntegralFar;
      break;
    default:
      break;
    }

    break;
  }
  case 2:
  {
    switch (_kIndexIncrement)
    {
    case 0:
      integralZ=fullIntegralNear;
      break;
    case 1:
      integralZ=fullIntegralFar;
      break;
    case 2:
      integralZ=0.0;
      break;
    case -1:
      integralZ=fullIntegralNear;
      break;
    case -2:
      integralZ=fullIntegralFar;
      break;
    default:
      break;
    }

    switch (_jIndexIncrement)
    {
    case 0:
      integralY=twiceHalfLowerNear;
      break;
    case 1:
      integralY=acrossNearFar;
      break;
    case 2:
      integralY=halfUpperIntegralFar;
      break;
    case -1:
      integralY=acrossNearFar;
      break;
    case -2:
      integralY=halfUpperIntegralFar;
      break;
    default:
      break;
    }

    switch (_iIndexIncrement)
    {
    case 0:
      integralX=twiceHalfLowerNear;
      break;
    case 1:
      integralX=acrossNearFar;
      break;
    case 2:
      integralX=halfUpperIntegralFar;
      break;
    case -1:
      integralX=acrossNearFar;
      break;
    case -2:
      integralX=halfUpperIntegralFar;
      break;
    default:
      break;
    }

    break;
  }
  default:
    break;
  }


  result=integralX*integralY*integralZ;


  return result;
}

//----------------------------------------------------------------------------------------------------------------------

float MathFunctions::calcTightQuadraticStencil(float _x)
{
  /*Tight Quadratic Stencil
  --------------------------------------------------------------------------------------------------
  N(x) = -x^2 +3/4                       if 0<=|x|<1/2
       = (1/2)*(x^2) - (3/2)|x| + 9/8      if 1/2<=|x|<3/2  NB! Paper says 1<=|x|<3/2 but assuming meant 1/2
       = 0                               otherwise
  --------------------------------------------------------------------------------------------------
  */

  float result=0.0;
  float absX=std::abs(_x);

  if (absX<0.5)
  {
    result=(-1.0)*(pow(_x,2));
    result+=(3.0/4.0);
  }

  else if (absX>=0.5 && absX<1.5)
  {
    result=0.5*(pow(_x,2));
//    result-=((3.0/2.0)*_x);
    result-=((3.0/2.0)*absX);
    result+=(9.0/8.0);
  }

  return result;
}

//----------------------------------------------------------------------------------------------------------------------

float MathFunctions::calcTightQuadraticStencil_Diff(float _x)
{
  /*Tight Quadratic Stencil
  --------------------------------------------------------------------------------------------------
  N(x) = -2x              if 0<=|x|<1/2
       = x - 3/2*signX    if 1/2<=|x|<3/2  NB! Paper says 1<=|x|<3/2 but assuming meant 1/2
       = 0                otherwise
  --------------------------------------------------------------------------------------------------
  */

  float result=0.0;
  float absX=std::abs(_x);

  float signX=MathFunctions::signFunction(_x);

  if (absX<0.5)
  {
    result=(-2.0*_x);
  }

  else if (absX>=0.5 && absX<1.5)
  {
    result=_x;
//    result-=(3.0/2.0);
    result-=((3.0/2.0)*signX);
  }

  return result;
}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::conjugateGradient(const Eigen::SparseMatrix<double> &_A, const Eigen::VectorXd &_B, Eigen::VectorXd &o_x, float _maxLoops, float _minResidual)
{
  /// @brief Function which uses Conjugate Gradient to solve Ax=b
  /// A has to be symmetric, definite and square.
  /// @todo Implement so can solve with preconditioner.

  //Initialise the Conjugate Gradient solver
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> conjGrad;

  //Insert A
  conjGrad.compute(_A);

  //Set max iterations and min residual
  conjGrad.setMaxIterations(_maxLoops);
  conjGrad.setTolerance(_minResidual);


  //Solve
  o_x=conjGrad.solve(_B);

  //Print out iteration number and error
  std::cout<<"Number of iterations: "<<conjGrad.iterations()<<"\n";
  std::cout<<"Error: "<<conjGrad.error()<<"\n";


}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::linearSystemSolve(const Eigen::Matrix3f &_A, const Eigen::Vector3f &_B, Eigen::Vector3f &o_x)
{
  /// @brief Function to solve linear system where there are no restrictions on A
  /// except that it is 3x3 as this method is quite slow.
  /// Uses ColPivHouseholderQR from Eigen library - Accuracy = 3/4 and speed 2/4.

  //Set up solver
  Eigen::ColPivHouseholderQR<Eigen::Matrix3f> qrSolver(_A);

  //Find solution x
  o_x=qrSolver.solve(_B);

}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::polarDecomposition(const Eigen::Matrix3f &_decomposeMatrix, Eigen::Matrix3f &o_R, Eigen::Matrix3f &o_S)
{
  /// @brief Calculates polar decomposition using singular value decomposition.
  /// decomposeMatrix=RS, decomposeMatrix=UXV*, R=UV*, S=VXV*
  /// Checks for determinant, as no polar decomposition when determinant is zero

  if (_decomposeMatrix.determinant()!=0.0)
  {
    //Set up matrices for the singular value decomposition
    Eigen::Matrix3f U;
    Eigen::Matrix3f X;
    Eigen::Matrix3f V;

    //Perform singular value decomposition
    singularValueDecomposition(_decomposeMatrix, U, X, V);

    //Calculate conjugate transpose of V; V*
    Eigen::Matrix3f V_conj;
    Eigen::Matrix3f V_conjTrans;
    V_conj=V.conjugate();
    V_conjTrans=V_conj.transpose();

    //Calculate R
    o_R=U*V_conjTrans;

    //Calculate S
    o_S=V*X*V_conjTrans;
  }
  else
  {
    throw std::invalid_argument("Determinant is zero");
  }

}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::singularValueDecomposition(const Eigen::Matrix3f &_decomposeMatrix, Eigen::Matrix3f &o_U, Eigen::Matrix3f &o_singularValues, Eigen::Matrix3f &o_V)
{
  /// @brief Uses Eigen library JacobiSVD. Set up so can only solve for 3x3 matrices.

  //Set up singular value decomposition solver
  Eigen::JacobiSVD<Eigen::Matrix3f> SVD_solver(_decomposeMatrix, Eigen::ComputeFullU | Eigen::ComputeFullV);

  //Store U and V
  o_U=SVD_solver.matrixU();
  o_V=SVD_solver.matrixV();

  //Calculate singular values
  Eigen::Vector3f singularValueVector=SVD_solver.singularValues();

  //Create diagonal matrix with singular values
  Eigen::Matrix3f singularValueMatrix;

  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      //Fill 3x3 matrix from 3x1 vector
      if (i==j)
      {
        singularValueMatrix(i,j)=singularValueVector(i);
      }
      else
      {
        singularValueMatrix(i,j)=0.0;
      }
    }
  }

  o_singularValues=singularValueMatrix;

}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::centralDifferenceGradient()
{

}

//----------------------------------------------------------------------------------------------------------------------

float MathFunctions::signFunction(float _x)
{
  float signValue=0.0;

  if (_x>0.0)
  {
    signValue=1.0;
  }

  else if (_x<0.0)
  {
    signValue=-1.0;
  }

  return signValue;
}

//----------------------------------------------------------------------------------------------------------------------

int MathFunctions::findMinVectorValue(const std::vector<int> &_vectorList)
{
  int vectorSize=_vectorList.size();

  int minValue=std::numeric_limits<int>::max();

  for (int i=0; i<vectorSize; i++)
  {
    int testValue=_vectorList.at(i);
    if (testValue<minValue && testValue>0)
    {
      minValue=testValue;
    }
  }

  return minValue;
}

//----------------------------------------------------------------------------------------------------------------------

Eigen::Matrix3f MathFunctions::matrixElementMultiplication(const Eigen::Matrix3f &_A, const Eigen::Matrix3f &_B)
{
  //Set return matrix
  Eigen::Matrix3f result;

  //Loop over all elements in matrix. Set up for 3x3 only, so 3^2 elements
  for (int i=0; i<9; i++)
  {
    result(i)=_A(i)*_B(i);
  }

  //Return result
  return result;
}

//----------------------------------------------------------------------------------------------------------------------
