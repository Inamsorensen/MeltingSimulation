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

void Grid::update(float _dt)
{

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcInterpolationWeights()
{

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::classifyCells()
{

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::setBoundaryVelocity()
{

}

//----------------------------------------------------------------------------------------------------------------------
