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

void Grid::setUpB_projectVelocity()
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
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::setUpA_projectVelocity()
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Loop over all cell centres

    Loop over all grid cell centres

        Check if face colliding and cell centre interior
          Get volume for each face
        Else use volume of cell

        Add, hence must be initialised to zero, GicGic' equal to (2/px +2/py + 2/pz) for i==j and -1/px/y/z otherwise x/y/z depends which value of i and j

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
