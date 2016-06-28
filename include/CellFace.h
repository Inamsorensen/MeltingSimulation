#ifndef CELLFACE
#define CELLFACE

#include <vector>

#include "Particle.h"

//------------------------------------------------------------------------------------------------------------------------------------------------------
/// @file CellFace.h
/// @brief Structure for a single cell face, ie in either x, y or z direction. Contains data required for the calculations, but need to be updated for
/// each time step as no data stored by grid. Also contains list of particles contained in cell to reduce number of loops in data transitions.
/// @author Ina M. Sorensen
/// @version 1.0
/// @date 27.06.16
///
/// @todo
//------------------------------------------------------------------------------------------------------------------------------------------------------


struct CellFace
{
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief List of particles inside the cell at current time step
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<Particle*> m_particles;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief List of interpolation weights for each of these particles.
  /// Used to transfer data from particles to grid
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<float> m_interpWeights_PartToGrid;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief List of interpolation weights differential for each of these particles.
  /// Used to transfer data from particles to grid
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<float> m_interpWeightsDiff_PartToGrid;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Collision state of cell. Is it colliding, empty or interior
  //----------------------------------------------------------------------------------------------------------------------
  State m_state;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Mass of cell
  //----------------------------------------------------------------------------------------------------------------------
  float m_mass;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Cell face velocity
  //----------------------------------------------------------------------------------------------------------------------
  float m_velocity;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Heat capacity for cell. Used in calculating temperature of cell
  //----------------------------------------------------------------------------------------------------------------------
  float m_heatCapacity;
};

#endif // CELLFACE

