#include "Grid.h"

#include <stdexcept>
#include <iostream>
#include <cmath>
#include <math.h>

//----------------------------------------------------------------------------------------------------------------------

Grid* Grid::m_instance=nullptr;

//----------------------------------------------------------------------------------------------------------------------

Grid::Grid(Eigen::Vector3f _origin, float _gridSize, int _noCells)
{
  /// @brief Sets grid variables defining size of grid and cells, and the grid origin.
  /// Uses this to set up the cell lists

  //Set up grid variables
  m_origin=_origin;
  m_gridSize=_gridSize;
  m_noCells=_noCells;
  m_cellSize=m_gridSize/((float)m_noCells);

//  m_cellCentres.reserve(m_noCells);
//  m_cellFacesX.reserve(m_noCells);
//  m_cellFacesY.reserve(m_noCells);
//  m_cellFacesZ.reserve(m_noCells);

  //Setup cell lists
  for (int i=0; i<pow(m_noCells,3); i++)
  {
    CellCentre* cellCentre=new CellCentre();
    CellFace* cellFaceX=new CellFace();
    CellFace* cellFaceY=new CellFace();
    CellFace* cellFaceZ=new CellFace();

    m_cellCentres.push_back(cellCentre);
    m_cellFacesX.push_back(cellFaceX);
    m_cellFacesY.push_back(cellFaceY);
    m_cellFacesZ.push_back(cellFaceZ);

  }


}

//----------------------------------------------------------------------------------------------------------------------

Grid::~Grid()
{
  /// @brief Delete lists of cell pointers and instance pointer

  int noCellCentresCurrent=m_cellCentres.size();
  int noCellFacesXCurrent=m_cellFacesX.size();
  int noCellFacesYCurrent=m_cellFacesY.size();
  int noCellFacesZCurrent=m_cellFacesZ.size();

  if (noCellCentresCurrent!=0 && noCellFacesXCurrent!=0 && noCellFacesYCurrent!=0 && noCellFacesZCurrent!=0)
  {
    if (noCellCentresCurrent==noCellFacesXCurrent==noCellFacesYCurrent==noCellFacesZCurrent)
    {
      for (int i=0; i<noCellCentresCurrent; i++)
      {
        delete m_cellCentres[i];
        delete m_cellFacesX[i];
        delete m_cellFacesY[i];
        delete m_cellFacesZ[i];
      }
    }
    else
    {
      for (int i=0; i<noCellCentresCurrent; i++)
      {
        delete m_cellCentres[i];
      }
      for (int i=0; i<noCellFacesXCurrent; i++)
      {
        delete m_cellFacesX[i];
      }
      for (int i=0; i<noCellFacesYCurrent; i++)
      {
        delete m_cellFacesY[i];
      }
      for (int i=0; i<noCellFacesZCurrent; i++)
      {
        delete m_cellFacesZ[i];
      }
    }
  }

  m_cellCentres.clear();
  m_cellFacesX.clear();
  m_cellFacesY.clear();
  m_cellFacesZ.clear();

  std::cout<<"Deleting grid\n";

}

//----------------------------------------------------------------------------------------------------------------------

Grid* Grid::createGrid(Eigen::Vector3f _origin, float _gridSize, int _noCells)
{
  /// @brief Creates grid from input variables. Only if grid has not already been created.

  if (m_instance==nullptr)
  {
    m_instance=new Grid(_origin, _gridSize, _noCells);
  }

  return m_instance;
}

//----------------------------------------------------------------------------------------------------------------------

Grid* Grid::getGrid()
{
  if (m_instance==nullptr)
  {
    throw std::invalid_argument("You need to create the grid first.");
  }

  return m_instance;

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::update(float _dt, bool _isFirstStep)
{
  /* Outline
  ---------------------------------------------------------------------------------------------------------------------
  Set m_dt=_dt

  Clear InterpolationData for each grid cell so all empty before start adding particles

  findParticleInCell - need to find out which particles are in which cells and their respective interp weight

  Get particle data to grid

  Classify cells

  Compute deviatoric force

  Set velocity due to collisions - Not sure why here

  Apply implicit velocity update

  Project velocity, ie. calc pressure

  Solve heat equation - do this in a separate thread?

  Update particle from grid

  ---------------------------------------------------------------------------------------------------------------
  */
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::clearCellData()
{
  /* Outline
  --------------------------------------------------------------------------------------------------------------
  Loop over all cells
    If m_interpolationData.size!=0, clear it
  --------------------------------------------------------------------------------------------------------------
  */
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::findParticleInCell()
{
  /* Outline
  ----------------------------------------------------------------------------------------------------
   Loop over all particles - This is where parallel should be inserted I think
   {
     Find position of particle in grid using getParticleGridCell
     This gives vector of i,j,k for cell

     Get neighbours i+-2, j+-2, k+-2. Ie loop over these

     Pass in cell i,j,k to calcInterpolationWeights
    }
  ------------------------------------------------------------------------------------------------------
  */
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcInterpolationWeights(Particle* _particle, int _i, int _j, int _k)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Calculate interpolation weights for the particle and the cell with given _i, _j,_k

  Calculate i*, ix, iy, iz for the cell in question

  Calculate x values for these

  Set up floats to recieve from maths functions

  Pass in x and get interpolation weights:
       cubicBspline
       cubicBspline_Diff
       cubicBspline_Integ
       tightQuadraticStencil
       tightQuadraticStencil_Diff

  Store interpolation weights in cell centre and face
  ------------------------------------------------------------------------------------------------------
  */

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::transferParticleData(bool _isFirstStep)
{
  /* Outline
  -----------------------------------------------------------------------------------------------------
  Loop over all cells - Parallelise here
  Could potentially set this so loop over all faces and centres, ie. use more threads

    Check that m_InterpolationData is not empty - If it is, set everything to zero?

    For each particle in the list calculate all variables
      m_i
      v_i
      kappa_i
      m_c
      J_c
      JE_c
      c
      T
      lambda^-1
      JP_c


     If first step in simulation, add to particle density

  If first step in simulation, loop over particles
    Calculate V_p


  -----------------------------------------------------------------------------------------------------
  */
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::classifyCells()
{
  /* Outline - Current setup might be time consuming since two loops. Could rectify this for bounding box, but not for level sets I think
  ----------------------------------------------------------------------------------------------------------------------
  Loop over all cell faces
    Check faces against collision
      For now use bounding box, so colliding if i<2||>n-2, j<2||>n-2, k<2||>n-2

  Loop over all cell centres
    Check if all 6 faces are colliding - If cell centre is i=n-1, j=n-1 or k=n-1 then only three faces
      If all colliding - colliding
      If not all and no particles - empty
      Otherwise - interior
    If colliding, also need to set temperature. Also need to set temp for empty


  ----------------------------------------------------------------------------------------------------------------------
  */
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::setBoundaryVelocity()
{
  /* Outline
  ---------------------------------------------------------------------------------------------------------------------
  Loop over all cells
    Set velocity to zero at colliding cell faces

  ---------------------------------------------------------------------------------------------------------------------
  */

}

//----------------------------------------------------------------------------------------------------------------------
