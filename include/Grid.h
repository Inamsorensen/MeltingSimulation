#ifndef GRID
#define GRID

#include <omp.h>

#include <vector>

#include <eigen3/Eigen/Core>

#include <ngl/Vec3.h>

#include "CellCentre.h"
#include "CellFace.h"
#include "Emitter.h"
#include "MathFunctions.h"

//------------------------------------------------------------------------------------------------------------------------------------------------------
/// @file Grid.h
/// @brief Grid class for the grid on which all calculations are done. Singleton class as only one grid for the calculation.
/// Contains list of pointers to all grid cell centres and faces.
/// @author Ina M. Sorensen
/// @version 1.0
/// @date 27.06.16
///
/// @done:Staggered grid
///
///
/// @todo Reserve memory for interpolation data vectors
///       Verify particle position in grid
///       Check that particles stored correctly in interp data
//------------------------------------------------------------------------------------------------------------------------------------------------------


class Grid
{
public:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Grid create instance, creates grid and returns instance
  /// @param [in] _origin is the position of the grid origin, set to lower, back left corner of grid
  /// @param [in] _gridSize is the size of one lenght of the grid. Grid is always cubic, so same lenght in all directions
  /// @param [in] _noCells is the number of grid cells in one direction. Same number in all directions
  //----------------------------------------------------------------------------------------------------------------------
  static Grid* createGrid(Eigen::Vector3f _originEdge, float _boundingBoxSize, int _noCells);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get instance of grid
  //----------------------------------------------------------------------------------------------------------------------
  static Grid* getGrid();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Destructor
  //----------------------------------------------------------------------------------------------------------------------
  ~Grid();

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get grid position. Different from bounding box as single layer of cells around bounding box for collision
  /// Returns position of lower back corner of grid, not the staggered position
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Vector3f getGridCornerPosition();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get cell size
  //----------------------------------------------------------------------------------------------------------------------
  float getGridCellSize(){return m_cellSize;}

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set surrounding temperatures; ambient temp and heat source temp.
  /// Reads in temperatures in celsius and sets them to kelvin for calculations
  //----------------------------------------------------------------------------------------------------------------------
  void setSurroundingTemperatures(float _ambientTemp, float _heatSourceTemp);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Find no of particles in each grid cell. Takes particle in emitter and checks positions against grid cells
  /// Returns vector containing the number of particles in each cell.
  /// @param [in] _emitter is used to get the list of particles to check for
  /// @param [out] o_listParticleNo: Stores the number of particles, one number for each cell. If _storeZero=true, then can
  /// access one grid cell using getVectorIndex from MathFunctions
  //----------------------------------------------------------------------------------------------------------------------
  void findNoParticlesInCells(Emitter* _emitter, std::vector<int> &o_listParticleNo);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Does all calculations for one time step. Updates velocity and temperature through force, pressure and
  /// temperature calculations
  //----------------------------------------------------------------------------------------------------------------------
  void update(float _dt, Emitter *_emitter, bool _isFirstStep, float _velocityContribAlpha, float _temperatureContribBeta);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get cell state for visualisation
  //----------------------------------------------------------------------------------------------------------------------
  inline State getCellState(int _cellIndex) const {return m_cellCentres[_cellIndex]->m_state;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get cell temperature for visualisation
  //----------------------------------------------------------------------------------------------------------------------
  inline float getCellTemperature(int _cellIndex) const {return m_cellCentres[_cellIndex]->m_temperature;}


private:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Constructor. Private for a singleton
  //----------------------------------------------------------------------------------------------------------------------
  Grid(Eigen::Vector3f _originEdge, float _boundingBoxSize, int _noCells);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Instance pointer
  //----------------------------------------------------------------------------------------------------------------------
  static Grid* m_instance;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Location of grid origin. Origin set to lower, back left corner.
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Vector3f m_origin;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Size of grid along one side. Grid always cubic so same length in all directions.
  //----------------------------------------------------------------------------------------------------------------------
  float m_gridSize;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Number of cells along one side. Same number along all directions
  //----------------------------------------------------------------------------------------------------------------------
  int m_noCells;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Total number of cells in grid.
  //----------------------------------------------------------------------------------------------------------------------
  int m_totNoCells;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Size of a single cell
  //----------------------------------------------------------------------------------------------------------------------
  float m_cellSize;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief List of cell centres. Each cell centre contains data for calculation
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<CellCentre*> m_cellCentres;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief List of cell faces, one in each directions of the face normals. Each cell face contains data for calculations.
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<CellFace*> m_cellFacesX;
  std::vector<CellFace*> m_cellFacesY;
  std::vector<CellFace*> m_cellFacesZ;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Number of empty cell centres and faces
  //----------------------------------------------------------------------------------------------------------------------
  int m_noCellCentres_Empty;
  int m_noCellFacesX_Empty;
  int m_noCellFacesY_Empty;
  int m_noCellFacesZ_Empty;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Number of interior cell centres and faces
  //----------------------------------------------------------------------------------------------------------------------
  int m_noCellCentres_Interior;
  int m_noCellFacesX_Interior;
  int m_noCellFacesY_Interior;
  int m_noCellFacesZ_Interior;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Simulation time step
  //----------------------------------------------------------------------------------------------------------------------
  float m_dt;


  //----------------------------------------------------------------------------------------------------------------------
  /// @brief External force on simulation. Set to gravity for now
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Vector3f m_externalForce;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief A matrix for deviatoric velocity calculation
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::MatrixXf m_Amatrix_deviatoric_X;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief B vector for deviatoric velocity calculation
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::VectorXf m_Bvector_deviatoric_X;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief A matrix for deviatoric velocity calculation
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::MatrixXf m_Amatrix_deviatoric_Y;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief B vector for deviatoric velocity calculation
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::VectorXf m_Bvector_deviatoric_Y;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief A matrix for deviatoric velocity calculation
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::MatrixXf m_Amatrix_deviatoric_Z;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief B vector for deviatoric velocity calculation
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::VectorXf m_Bvector_deviatoric_Z;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Ambient temperature; temperature of the surrounding air. In Kelvin
  //----------------------------------------------------------------------------------------------------------------------
  float m_ambientTemperature;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Heat source temperature. In Kelvin
  //----------------------------------------------------------------------------------------------------------------------
  float m_heatSourceTemperature;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Threshold for number of particles that must be affecting cell. Otherwise get too small mass.
  //----------------------------------------------------------------------------------------------------------------------
  int m_noParticlesThreshold;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Determine whether should use implicit or explicit intergration for deviatoric velocity
  //----------------------------------------------------------------------------------------------------------------------
  float m_isImplictIntegration;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Clear list of InterpolationData
  //----------------------------------------------------------------------------------------------------------------------
  void clearCellData();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Find which cells have particles in or near them such that interpolation weight will not be zero
  //----------------------------------------------------------------------------------------------------------------------
  void findParticleContributionToCell(Emitter *_emitter);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate interpolation weights for transitions Particle-Grid and Grid-Particle
  //----------------------------------------------------------------------------------------------------------------------
  void calcInterpolationWeights(Particle *_particle, int _i, int _j, int _k);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Transfer particle data to grid
  //----------------------------------------------------------------------------------------------------------------------
  void transferParticleData(Emitter *_emitter);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate initial particle volumes
  //----------------------------------------------------------------------------------------------------------------------
  void calcInitialParticleVolumes(Emitter* _emitter);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Verify whether cell centres and faces are colliding, empty or interior
  /// @todo Change to switch/case statements instead of if statements
  //----------------------------------------------------------------------------------------------------------------------
  void classifyCells();

  //NEW INTERPOLATION AND DEVIATORIC CALC SETUP - 14.08.16
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Interpolate data from particles to grid
  /// @brief Calculates variables for grid
  /// @brief Calculates deviatoric force, B component and A components for deviatoric velocity calculations
  //----------------------------------------------------------------------------------------------------------------------
  void interpolateParticleToGrid(Emitter* _emitter, bool _isFirstStep);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate interpolation weights based on cubic B Spline for given particle and cell
  //----------------------------------------------------------------------------------------------------------------------
  void calcWeight_cubicBSpline(Eigen::Vector3f _particlePosition, int _iIndex, int _jIndex, int _kIndex,
                                    float &o_weightCentre, float &o_weightFaceX, float &o_weightFaceY, float &o_weightFaceZ);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate interpolation weights differentiated based on cubic B Spline for given particle and cell
  //----------------------------------------------------------------------------------------------------------------------
  void calcWeight_cubicBSpline_Diff(Eigen::Vector3f _particlePosition, int _iIndex, int _jIndex, int _kIndex,
                                    Eigen::Vector3f &o_weightFaceX, Eigen::Vector3f &o_weightFaceY, Eigen::Vector3f &o_weightFaceZ);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Verify whether cell centres and faces are colliding, empty or interior
  /// @todo Change to switch/case statements instead of if statements
  //----------------------------------------------------------------------------------------------------------------------
  void classifyCells_New();

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate deviatoric force, B component and A row component contributions from a particle
  //----------------------------------------------------------------------------------------------------------------------
  void calcDeviatoricContributions(Particle *_particle, int _cellIndex, int _iIndex, int _jIndex, int _kIndex, float _weightX, float _weightY, float _weightZ);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate force due to deviatoric stress
  //----------------------------------------------------------------------------------------------------------------------
  float calcDeviatoricForce_New(Particle *_particle, Eigen::Vector3f _eVector, Eigen::Vector3f _weightDiff);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate B in Ax=B for the implicit calculation of velocity.
  /// b_i=v_i + (dt/m_i)*f_i + dt*g*eVector*sumWeight
  //----------------------------------------------------------------------------------------------------------------------
  float calcBComponent_DeviatoricVelocity_New(Particle *_particle, Eigen::Vector3f _eVector, float _weight, float _deviatoricForce, float _massFace);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get A components of one row for deviatoric velocity calculation
  //----------------------------------------------------------------------------------------------------------------------
  void calcAComponent_DeviatoricVelocity_New(Particle* _particle, int _cellIndex_column, Eigen::Vector3f _weightDiff_FaceX_column, Eigen::Vector3f _weightDiff_FaceY_column, Eigen::Vector3f _weightDiff_FaceZ_column);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculates Ap matrix which is used to calculate the A matrix components for implicit deviatoric velocity integration
  /// @brief Ap=d2Y_hat/dFE2 : eVector*weight_diff_trans*deformGradElastic
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Matrix3f calcApComponent_DeviatoricVelocity_New(Particle* _particle, Eigen::Vector3f _weight_diff_column, Eigen::Vector3f _eVector);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate deltaR
  /// @brief This is the delta of the rotation matrix R in polar decomposition of the deformation gradient: F=RS
  /// @brief where F is defined as J^(-1/dimension)*F, also called deviatoric F.
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Matrix3f calculate_dR_New(const Eigen::Matrix3f &_deltaDeformElastic_Deviatoric, const Eigen::Matrix3f &_R_deformElastic_Deviatoric, const Eigen::Matrix3f &_S_deformElastic_Deviatoric);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Explicitly update velocity. v_new=bcomponent
  //----------------------------------------------------------------------------------------------------------------------
  void explicitUpdate_DeviatoricVelocity_New();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Implicitly update velocity.
  //----------------------------------------------------------------------------------------------------------------------
  void implicitUpdate_DeviatoricVelocity_New();

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Update particle data from grid
  //----------------------------------------------------------------------------------------------------------------------
  void updateParticleFromGrid_New(Emitter *_emitter, float _velocityContribAlpha, float _tempContribBeta);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate interpolation weights based on cubic B Spline for given particle and cell
  //----------------------------------------------------------------------------------------------------------------------
  void calcWeight_tightQuadraticStencil(Eigen::Vector3f _particlePosition, int _iIndex, int _jIndex, int _kIndex,
                                    float &o_weightCentre, float &o_weightFaceX, float &o_weightFaceY, float &o_weightFaceZ);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate interpolation weights differentiated based on cubic B Spline for given particle and cell
  //----------------------------------------------------------------------------------------------------------------------
  void calcWeight_tightQuadraticStencil_Diff(Eigen::Vector3f _particlePosition, int _iIndex, int _jIndex, int _kIndex,
                                             Eigen::Vector3f &o_weightFaceX, Eigen::Vector3f &o_weightFaceY, Eigen::Vector3f &o_weightFaceZ);

  //END NEW INTERPOLATION AND DEVIATORIC CALC SETUP - 14.08.16

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate force due to deviatoric stress
  //----------------------------------------------------------------------------------------------------------------------
  float calcDeviatoricForce(Particle *_particle, Eigen::Vector3f _eVector, Eigen::Vector3f _weightDiff);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate velocity after external forces and deviatoric stress has been applied
  /// @todo Set boundary velocities
  //----------------------------------------------------------------------------------------------------------------------
  void calcDeviatoricVelocity();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate B in Ax=B for the implicit calculation of velocity.
  /// b_i=v_i + (dt/m_i)*f_i + dt*g*eVector*sumWeight
  //----------------------------------------------------------------------------------------------------------------------
  float calcBComponent_DeviatoricVelocity(CellFace *_cellFace, Eigen::Vector3f _eVector);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get A components of one row for deviatoric velocity calculation
  //----------------------------------------------------------------------------------------------------------------------
  void calcAComponent_DeviatoricVelocity(int _cellIndex, int _noParticlesFaceX, int _noParticlesFaceY, int _noParticlesFaceZ, Eigen::MatrixXf &o_AX, Eigen::MatrixXf &o_AY, Eigen::MatrixXf &o_AZ);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculates component of A matrix for Ax=b. In this case have (I+A)x=b where I will not be included in the A component
  //----------------------------------------------------------------------------------------------------------------------
  float calcAValue_DeviatoricVelocity(Particle* _particle, Eigen::Vector3f _weight_i_diff, Eigen::Vector3f _weight_j_diff, Eigen::Vector3f _eVector);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Explicitly update velocity. v_new=bcomponent
  //----------------------------------------------------------------------------------------------------------------------
  void explicitUpdateVelocity(int _cellIndex, float _velocityX, float _velocityY, float _velocityZ);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Implicitly update velocity.
  //----------------------------------------------------------------------------------------------------------------------
  void implicitUpdateVelocity(const Eigen::MatrixXf &_A_X, const Eigen::VectorXf &_bVector_X, const Eigen::MatrixXf &_A_Y, const Eigen::VectorXf &_bVector_Y, const Eigen::MatrixXf &_A_Z, const Eigen::VectorXf &_bVector_Z);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate rotation matrix R in polar decomposition of the deformation gradient: F=RS
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Matrix3f calculate_dR(const Eigen::Matrix3f &_deltaDeformElastic_Deviatoric, const Eigen::Matrix3f &_R_deformElastic_Deviatoric, const Eigen::Matrix3f &_S_deformElastic_Deviatoric);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Search list of particles to find same particle in two cells
  //----------------------------------------------------------------------------------------------------------------------
  void searchCellsForCommonParticle(unsigned int _particleId, CellFace* _cellFace, unsigned int &o_particleIndexInFace, bool &o_isFound);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set boundary velocity. Set as stick on collision, ie. zero velocity for colliding faces
  //----------------------------------------------------------------------------------------------------------------------
  void setBoundaryVelocity();

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Project velocity to find velocity which agrees with pressure calculations such that incompressibility can
  /// be accounted for
  ///
  /// @todo Read through properly
  ///       Check ways of parallelising
  //----------------------------------------------------------------------------------------------------------------------
  void projectVelocity();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set up B in Ax=B for Poisson equation which solves for pressure
  //----------------------------------------------------------------------------------------------------------------------
  float calcBComponent_projectVelocity(int _cellIndex, int _iIndex, int _jIndex, int _kIndex);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set up A in Ax=B for Poisson equation which solves for pressure
  //----------------------------------------------------------------------------------------------------------------------
  void calcAComponent_projectVelocity(int _cellIndex, int _iIndex, int _jIndex, int _kIndex, Eigen::SparseMatrix<double> &o_A, Eigen::MatrixXf &o_A_test);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate exact volume of cell at boundaries. If not done, then this volume will be too small, and lead to
  /// errors
  //----------------------------------------------------------------------------------------------------------------------
  void calcFaceDensities(int _cellIndex);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate cell face control volume
  //----------------------------------------------------------------------------------------------------------------------
  float calcFaceVolume(int _iIndex, int _jIndex, int _kIndex, int _cellFaceDirection);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate new temperature values
  //----------------------------------------------------------------------------------------------------------------------
  void calcTemperature();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set up B in Ax=B to solve for temperature
  //----------------------------------------------------------------------------------------------------------------------
  float calcBComponent_temperature(int _cellIndex);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set up A in Ax=B to solve for temperature
  //----------------------------------------------------------------------------------------------------------------------
  void calcAComponent_temperature(int _cellIndex, int _iIndex, int _jIndex, int _kIndex, Eigen::SparseMatrix<double> &o_A);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Update particle data from grid
  //----------------------------------------------------------------------------------------------------------------------
  void updateParticleFromGrid(float _velocityContribAlpha, float _tempContribBeta);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Update particle position directly
  //----------------------------------------------------------------------------------------------------------------------
  void updateParticlePositionDirectly(float _velocityContribAlpha, int _cellIndex);


};

#endif // GRID

