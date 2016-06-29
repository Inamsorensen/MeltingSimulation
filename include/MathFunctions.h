#ifndef MATHFUNCTIONS
#define MATHFUNCTIONS

#include <vector>

#include <ngl/Vec3.h>
#include <ngl/Mat3.h>

//------------------------------------------------------------------------------------------------------------------------------------------------------
/// @file MathFunctions.h
/// @brief Structure which contains general mathematical calculations. Most of these will be using the Eigen library.
/// @author Ina M. Sorensen
/// @version 1.0
/// @date 27.06.16
///
/// @todo
//------------------------------------------------------------------------------------------------------------------------------------------------------


struct MathFunctions
{
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get vector index from (i,j,k) cell index
  /// @param [in] Cell index in 3d (i,j,k)
  /// @param [in] _noCells is the total number of cells in grid
  //----------------------------------------------------------------------------------------------------------------------
  int getVectorIndex(int i, int j, int k, int _noCells);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get cell index of a particle
  /// @param [in] Position of particle
  //----------------------------------------------------------------------------------------------------------------------
  ngl::Vec3 getParticleGridCell(ngl::Vec3 _particlePosition, float _cellSize, ngl::Vec3 _gridOrigin);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Returns value of cubic B-spline
  /// @param [in] Position in a single direction, x
  //----------------------------------------------------------------------------------------------------------------------
  float calcCubicBSpline(float _x);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Returns value of differentiated cubic B-spline
  /// @param [in] Position in a single direction, x
  //----------------------------------------------------------------------------------------------------------------------
  float calcCubicBSpline_Diff(float _x);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Returns value of integrated cubic B-spline
  /// @param [in] Position in a single direction, x
  //----------------------------------------------------------------------------------------------------------------------
  float calcCubicBSpline_Integ(float _x);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Returns value of tight quadratic stencil
  /// @param [in] Position in a single direction, x
  //----------------------------------------------------------------------------------------------------------------------
  float calcTightQuadraticStencil(float _x);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Returns value of differentiated tight quadratic stencil
  /// @param [in] Position in a single direction, x
  //----------------------------------------------------------------------------------------------------------------------
  float calcTightQuadraticStencil_Diff(float _x);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Solves Ax=B using conjugate residual method
  /// @param [in] _A and _B which are a 2 and 1 dimensional matrix respectively.
  /// @param [in] _maxLoops is the max number of loops the method will do unless _minResidual is met first.
  /// @param [in] _x0 is a 1 dimensional vector giving the first guess at the solution
  /// @param[out] o_x is the solution
  //----------------------------------------------------------------------------------------------------------------------
  void conjugateResidual(std::vector<float> _A, std::vector<float> _B, std::vector<float> _x0, std::vector<float> o_x, float _maxLoops, float _minResidual);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Solves Ax=B using conjugate gradient method
  /// @param [in] _A and _B which are a 2 and 1 dimensional matrix respectively.
  /// @param [in] _maxLoops is the max number of loops the method will do unless _minResidual is met first.
  /// @param [in] _x0 is a 1 dimensional vector giving the first guess at the solution
  /// @param[out] o_x is the solution
  //----------------------------------------------------------------------------------------------------------------------
  void conjugateGradient(std::vector<float> _A, std::vector<float> _B, std::vector<float> _x0, std::vector<float> o_x, float _maxLoops, float _minResidual);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Solves Ax=B for all possible matrices A. Will use a method in Eigen that is slow, so only used for small matrices
  /// @param [in] _A and _B which are a 2 and 1 dimensional matrix respectively.
  /// @param[out] o_x is the solution
  //----------------------------------------------------------------------------------------------------------------------
  void linearSystemSolve(std::vector<float> _A, std::vector<float> _B, std::vector<float> o_x);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Polar decomposition
  /// @param [in] _decomposeMatrix is the matrix to be polar decomposed
  /// @param [out] o_R is the rotation matrix, o_S is the stretch matrix
  //----------------------------------------------------------------------------------------------------------------------
  void polarDecomposition(ngl::Mat3 _decomposeMatrix, ngl::Mat3 o_R, ngl::Mat3 o_S);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Singular value decomposition
  /// @param [in] _decomposeMatrix is the matrix to be decomposed
  /// @param [out] o_U and o_V are the left and right singular vectors respectively
  /// @param [out] o_singularValues is a diagonal matrix containing the singular values.
  //----------------------------------------------------------------------------------------------------------------------
  void singularValueDecomposition(ngl::Mat3 _decomposeMatrix, ngl::Mat3 o_U, ngl::Mat3 o_singularValues, ngl::Mat3 o_V);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculates the central difference gradient
  //----------------------------------------------------------------------------------------------------------------------
  void centralDifferenceGradient();

};

#endif // MATHFUNCTIONS

