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
  /// @brief Simulation time step
  //----------------------------------------------------------------------------------------------------------------------
  float m_dt;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief External force on simulation. Set to gravity for now
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Vector3f m_externalForce;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Ambient temperature; temperature of the surrounding air. In Kelvin
  //----------------------------------------------------------------------------------------------------------------------
  float m_ambientTemperature;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Heat source temperature. In Kelvin
  //----------------------------------------------------------------------------------------------------------------------
  float m_heatSourceTemperature;

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
  float calcBComponent_DeviatoricVelocity(float _velocity, float _mass, float _deviatoricForce, float _sumWeight, Eigen::Vector3f _eVector);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculates component of A matrix for Ax=b. In this case have (I+A)x=b where I will not be included in the A component
  //----------------------------------------------------------------------------------------------------------------------
  float calcAComponent_DeviatoricVelocity(Particle* _particle, Eigen::Vector3f _weight_i_diff, Eigen::Vector3f _weight_j_diff, Eigen::Vector3f _eVector, float _mass_i);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Explicitly update velocity. v_new=bcomponent
  //----------------------------------------------------------------------------------------------------------------------
  void explicitUpdateVelocity(int _cellIndex, float _velocityX, float _velocityY, float _velocityZ);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Implicitly update velocity.
  //----------------------------------------------------------------------------------------------------------------------
  void implicitUpdateVelocity(const Eigen::VectorXf &_bVector_X, const Eigen::VectorXf &_bVector_Y, const Eigen::VectorXf &_bVector_Z);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate rotation matrix R in polar decomposition of the deformation gradient: F=RS
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Matrix3f calculate_dR(const Eigen::Matrix3f &_deltaDeformElastic_Deviatoric, const Eigen::Matrix3f &_R_deformElastic_Deviatoric, const Eigen::Matrix3f &_S_deformElastic_Deviatoric);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set boundary velocity. Set as stick on collision, ie. zero velocity for colliding faces
  //----------------------------------------------------------------------------------------------------------------------
  void setBoundaryVelocity();

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Project velocity to find velocity which agrees with pressure calculations such that incompressibility can
  /// be accounted for
  //----------------------------------------------------------------------------------------------------------------------
  void projectVelocity();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set up B in Ax=B for Poisson equation which solves for pressure
  //----------------------------------------------------------------------------------------------------------------------
  void setUpB_projectVelocity();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set up A in Ax=B for Poisson equation which solves for pressure
  //----------------------------------------------------------------------------------------------------------------------
  void setUpA_projectVelocity();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set boundary pressure values. Set using Neumann boundary condition.
  //----------------------------------------------------------------------------------------------------------------------
  void setBoundaryPressure();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate exact volume of cell at boundaries. If not done, then this volume will be too small, and lead to
  /// errors
  //----------------------------------------------------------------------------------------------------------------------
  void calcBoundaryVolume();

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate new temperature values
  //----------------------------------------------------------------------------------------------------------------------
  void calcTemperature();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set up B in Ax=B to solve for temperature
  //----------------------------------------------------------------------------------------------------------------------
  void setUpB_temperature();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set up A in Ax=B to solve for temperature
  //----------------------------------------------------------------------------------------------------------------------
  void setUpA_temperature();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set boundary values for temperature. Set using Neumann boundary condition.
  //----------------------------------------------------------------------------------------------------------------------
  void setBoundaryTemperature();

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Update particle data from grid
  //----------------------------------------------------------------------------------------------------------------------
  void updateParticleFromGrid(float _velocityContribAlpha, float _tempContribBeta);


};

#endif // GRID

