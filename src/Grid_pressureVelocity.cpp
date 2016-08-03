#include "Grid.h"

//----------------------------------------------------------------------------------------------------------------------

void Grid::projectVelocity()
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Set up eigen matrices to store values in

  Set up B

  Set up A

  set Boundary conditions

  Solve using conjugate gradient

  Set result to cells

  ----------------------------------------------------------------------------------------------------------------
  */
}

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcBComponent_projectVelocity(int _cellIndex)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Loop over all grid cells

    For empty cells set equal to zero? Or should this be done afterwards so cells here definitely zero?

    Calculate central difference of velocity for each face and add together

    Calculate constant

    Add and insert into B vector


  ----------------------------------------------------------------------------------------------------------------
  */

    //Calculate constant
    float detDeformationGradElastic=m_cellCentres[_cellIndex]->m_detDeformationGradElastic;
    float constant=(detDeformationGradElastic-1.0);
    constant*=(1.0/m_dt*detDeformationGradElastic);
    constant*=(-1.0);

    //Get indices of current cell
    int iIndex=m_cellCentres[_cellIndex]->m_iIndex;
    int jIndex=m_cellCentres[_cellIndex]->m_jIndex;
    int kIndex=m_cellCentres[_cellIndex]->m_kIndex;

    //Get index of i+1,j,k and i,j+1,k and i,j,k+1
    int index_i1jk=MathFunctions::getVectorIndex(iIndex+1, jIndex, kIndex, m_noCells);
    int index_ij1k=MathFunctions::getVectorIndex(iIndex, jIndex+1, kIndex, m_noCells);
    int index_ijk1=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex+1, m_noCells);

    //Get face velocities for all faces surrounding cell centre
    float velocityX_forward=m_cellFacesX[index_i1jk]->m_velocity;
    float velocityX_backward=m_cellFacesX[_cellIndex]->m_velocity;
    float velocityY_forward=m_cellFacesY[index_ij1k]->m_velocity;
    float velocityY_backward=m_cellFacesY[_cellIndex]->m_velocity;
    float velocityZ_forward=m_cellFacesZ[index_ijk1]->m_velocity;
    float velocityZ_backward=m_cellFacesZ[_cellIndex]->m_velocity;

    //Calculate central gradient stencil
    float centralGradient=(1.0/m_cellSize)*((velocityX_forward-velocityX_backward)+(velocityY_forward-velocityY_backward)+(velocityZ_forward-velocityZ_backward));

    //Calculate B component
    float BComponent=constant-centralGradient;

    return BComponent;

}

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcAComponent_projectVelocity(int _cellIndex)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Loop over all cell centres

    Loop over all grid cell centres

        Check if face colliding and cell centre interior
          Get volume for each face
        Else use volume of cell

        Add, hence must be initialised to zero, GicGic' equal to (2/px +2/py + 2/pz) for i==j and -1/px/y/z otherwise x/y/z depends which value of i and j
        Need to set to zero for certain j values since only nearest neighbours have effect

        Check if colliding as need to take 1/px/y/z away from ii and set ij to zero

     End loop

     Multiply by dt/density_{i}

    End loop

  End loop

  Calculate constant

  Add constant*I - Not sure if this is faster than another loop

  ----------------------------------------------------------------------------------------------------------------
  */
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::setBoundaryPressure()
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Loop over all grid cells
    If empty
      Set pressure to zero

  ----------------------------------------------------------------------------------------------------------------
  */
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcBoundaryVolume()
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Leave for now I think

  ----------------------------------------------------------------------------------------------------------------
  */
}
