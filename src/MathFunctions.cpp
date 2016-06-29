#include "MathFunctions.h"

#include <stdexcept>
#include <cmath>
#include <math.h>

//----------------------------------------------------------------------------------------------------------------------

int MathFunctions::getVectorIndex(int i, int j, int k, int _noCells)
{
  int vectorIndex=i+(_noCells*j)+(_noCells*_noCells*k);
  return vectorIndex;
}

//----------------------------------------------------------------------------------------------------------------------

ngl::Vec3 MathFunctions::getParticleGridCell(ngl::Vec3 _particlePosition, float _cellSize, ngl::Vec3 _gridOrigin)
{
  ngl::Vec3 index;

  //Find grid indices from particle position.
  ngl::Vec3 indexParticle=(1.0/_cellSize)*(_particlePosition-_gridOrigin);

  //Find which cell the particle is in
  index.m_x=floor(indexParticle.m_x);
  index.m_y=floor(indexParticle.m_y);
  index.m_z=floor(indexParticle.m_z);

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

void MathFunctions::conjugateResidual(std::vector<float> _A, std::vector<float> _B, std::vector<float> _x0, std::vector<float> o_x, float _maxLoops, float _minResidual)
{

}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::conjugateGradient(std::vector<float> _A, std::vector<float> _B, std::vector<float> _x0, std::vector<float> o_x, float _maxLoops, float _minResidual)
{

}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::linearSystemSolve(std::vector<float> _A, std::vector<float> _B, std::vector<float> o_x)
{

}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::polarDecomposition(ngl::Mat3 _decomposeMatrix, ngl::Mat3 o_R, ngl::Mat3 o_S)
{

}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::singularValueDecomposition(ngl::Mat3 _decomposeMatrix, ngl::Mat3 o_U, ngl::Mat3 o_singularValues, ngl::Mat3 o_V)
{

}

//----------------------------------------------------------------------------------------------------------------------

void MathFunctions::centralDifferenceGradient()
{

}

//----------------------------------------------------------------------------------------------------------------------
