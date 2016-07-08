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

float MathFunctions::calcCubicBSpline_Integ(float _x)
{
  float result=0.0;
  return result;
}

//----------------------------------------------------------------------------------------------------------------------

float MathFunctions::calcTightQuadraticStencil(float _x)
{
  /*Tight Quadratic Stencil
  --------------------------------------------------------------------------------------------------
  N(x) = -x^2 +3/4                       if 0<=|x|<1/2
       = (1/2)*(x^2) - (3/2)x + 9/8      if 1/2<=|x|<3/2  NB! Paper says 1<=|x|<3/2 but assuming meant 1/2
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
    result-=((3.0/2.0)*_x);
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
       = x - 3/2          if 1/2<=|x|<3/2  NB! Paper says 1<=|x|<3/2 but assuming meant 1/2
       = 0                otherwise
  --------------------------------------------------------------------------------------------------
  */

  float result=0.0;
  float absX=std::abs(_x);

  if (absX<0.5)
  {
    result=(-2.0*_x);
  }

  else if (absX>=0.5 && absX<1.5)
  {
    result=_x;
    result-=(3.0/2.0);
  }

  return result;
}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::conjugateGradient(std::vector<float> *_A, std::vector<float> *_B, std::vector<float> *_x0, std::vector<float> *o_x, float _maxLoops, float _minResidual)
{
  /// @brief Function which uses Conjugate Gradient to solve Ax=b
  /// A has to be symmetric, definite and square.
  /// @todo Implement so can solve with guess.

  //Find size of the matrices
  int noRows=_B->size();

  //In case A is not a square matrix
  int noColumns=_A->size()/noRows;

  //Verify that _A is a square matrix, if not add line(s) of zero to A and b
  if (noRows!=noColumns)
  {
    throw std::invalid_argument("Only works for square matrix A");
  }

  //Set up matrices to be used by Eigen function
  Eigen::VectorXd b(noRows);
  Eigen::VectorXd x(noColumns);
  Eigen::SparseMatrix<double> A(noRows, noColumns);


  //Fill A and b
  for (int i=0; i<noRows; i++)
  {
    b(i)=_B->at(i);

    for (int j=0; j<noColumns; j++)
    {
      //Get vector index from i and j
      int index=j+(i*noColumns);

      A.insert(i,j)=_A->at(index);
    }
  }

  //Initialise the Conjugate Gradient solver
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> conjGrad;

  //Insert A
  conjGrad.compute(A);

  //Set max iterations and min residual
  conjGrad.setMaxIterations(_maxLoops);
  conjGrad.setTolerance(_minResidual);


  //Solve
  x=conjGrad.solve(b);

  //Print out iteration number and error
  std::cout<<"Number of iterations: "<<conjGrad.iterations()<<"\n";
  std::cout<<"Error: "<<conjGrad.error()<<"\n";

  //Put solution x into o_x for the return
  for (int i=0; i<noColumns; i++)
  {
    float value=x(i);
    o_x->push_back(value);
  }

}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::linearSystemSolve(std::vector<float>* _A, std::vector<float>* _B, std::vector<float>* o_x)
{
  /// @brief Function to solve linear system where there are no restrictions on A
  /// except that it is 3x3 as this method is quite slow.
  /// Uses ColPivHouseholderQR from Eigen library - Accuracy = 3/4 and speed 2/4.

  //Set up Eigen matrices
  Eigen::Matrix3f A;
  Eigen::Vector3f b;

  //Check that 3x3 matrix
  int sizeA=_A->size();
  int noRows=_B->size();
  int noColumns=sizeA/noRows;

  if (noRows!=3 || noColumns!=3)
  {
    throw std::invalid_argument("A must be 3x3 and b 3x1 matrices");
  }

  //Read in values
  for (int i=0; i<3; i++)
  {
    //Fill b
    b(i)=_B->at(i);

    for (int j=0; j<3; j++)
    {
      //Find index for _A vector
      int index=j+(i*3);

      //Fill A
      A(i,j)=_A->at(index);
    }
  }

  //Set up solver
  Eigen::ColPivHouseholderQR<Eigen::Matrix3f> qrSolver(A);

  //Find solution x
  Eigen::Vector3f x;
  x=qrSolver.solve(b);


  //Read solution into o_x
  for (int i=0; i<3; i++)
  {
    o_x->push_back(x[i]);
  }

}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::polarDecomposition(Eigen::Matrix3f *_decomposeMatrix, Eigen::Matrix3f *o_R, Eigen::Matrix3f *o_S)
{
  /// @brief Calculates polar decomposition using singular value decomposition.
  /// decomposeMatrix=RS, decomposeMatrix=UXV*, R=UV*, S=VXV*
  /// Checks for determinant, as no polar decomposition when determinant is zero

  if (_decomposeMatrix->determinant()!=0.0)
  {
    //Set up matrices for the singular value decomposition
    Eigen::Matrix3f U;
    Eigen::Matrix3f X;
    Eigen::Matrix3f V;

    //Perform singular value decomposition
    singularValueDecomposition(_decomposeMatrix, &U, &X, &V);

    //Calculate conjugate transpose of V; V*
    Eigen::Matrix3f V_conj;
    Eigen::Matrix3f V_conjTrans;
    V_conj=V.conjugate();
    V_conjTrans=V_conj.transpose();

    //Calculate R
    *o_R=U*V_conjTrans;

    //Calculate S
    *o_S=V*X*V_conjTrans;
  }
  else
  {
    throw std::invalid_argument("Determinant is zero");
  }

}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::singularValueDecomposition(Eigen::Matrix3f *_decomposeMatrix, Eigen::Matrix3f *o_U, Eigen::Matrix3f *o_singularValues, Eigen::Matrix3f *o_V)
{
  /// @brief Uses Eigen library JacobiSVD. Set up so can only solve for 3x3 matrices.

  //Set up singular value decomposition solver
  Eigen::JacobiSVD<Eigen::Matrix3f> SVD_solver(*_decomposeMatrix, Eigen::ComputeFullU | Eigen::ComputeFullV);

  //Store U and V
  *o_U=SVD_solver.matrixU();
  *o_V=SVD_solver.matrixV();

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

  *o_singularValues=singularValueMatrix;

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

int MathFunctions::findMinVectorValue(std::vector<int> *_vectorList)
{
  int vectorSize=_vectorList->size();

  int minValue=std::numeric_limits<int>::max();

  for (int i=0; i<vectorSize; i++)
  {
    int testValue=_vectorList->at(i);
    if (testValue<minValue && testValue>0)
    {
      minValue=testValue;
    }
  }

  return minValue;
}

//----------------------------------------------------------------------------------------------------------------------
