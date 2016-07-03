#include "MathFunctions.h"

#include <stdexcept>
#include <cmath>
#include <math.h>
#include <iostream>

//----------------------------------------------------------------------------------------------------------------------

int MathFunctions::getVectorIndex(int i, int j, int k, int _noCells)
{
  int vectorIndex=i+(_noCells*j)+(_noCells*_noCells*k);
  return vectorIndex;
}

//----------------------------------------------------------------------------------------------------------------------

Eigen::Vector3f MathFunctions::getParticleGridCell(Eigen::Vector3f _particlePosition, float _cellSize, Eigen::Vector3f _gridOrigin)
{
  Eigen::Vector3f index;

  //Find grid indices from particle position.
  Eigen::Vector3f indexParticle=(1.0/_cellSize)*(_particlePosition-_gridOrigin);

  //Find which cell the particle is in
  index(0)=floor(indexParticle(0));
  index(1)=floor(indexParticle(1));
  index(2)=floor(indexParticle(2));

  return index;
}

//----------------------------------------------------------------------------------------------------------------------

float MathFunctions::calcCubicBSpline(float _x)
{
  float result=0.0;
  return result;
}

//----------------------------------------------------------------------------------------------------------------------

float MathFunctions::calcCubicBSpline_Diff(float _x)
{
  float result=0.0;
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
  float result=0.0;
  return result;
}

//----------------------------------------------------------------------------------------------------------------------

float MathFunctions::calcTightQuadraticStencil_Diff(float _x)
{
  float result=0.0;
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
