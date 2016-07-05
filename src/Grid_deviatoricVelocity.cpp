#include "Grid.h"

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcDeviatoricForce()
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  f_{i}=-sum_p(V_{p}*e_{a(i)}^T*dYdFE_{p}:B_{kmij}*FE_{p}^T*cubicBSpline_Diff_{ip}

  dYdFE_{p}=2*lame_mu_{p}*JE^a*FE_{p} - 2*lame_mu_{p}*RE_{p} where a=-1/d where d=3? for dimensions

  B_{kmij}=JE^a*I + a*JE^a*FE_{p}^-T*FE_{p}

  Loop over faces - Could thread here and split further for next loop

  Loop over grid cells - Parallelise here
  Think calculate force for each of the faces but the only thing that changes is e_{a(i)}^T and cubicBSpline_Diff
  Depends on whether its the same particles for all faces

    Loop over particles in interpolationData

      calculate polar decomposition -  Actually should probably be done for each particle before grid as will need this and singular values later

      calc B_{kmij}

      calc dYdFE_{p}

      special multiply dYdFE_{p} and B_{kmij}

      multiply with V_{p}

      Loop over faces
        Multiply everything with weight and e specific to face

        Add to face's force

  ----------------------------------------------------------------------------------------------------------------
  */
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcDeviatoricVelocity()
{

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::setUpB_DeviatoricVelocity()
{
  /* Outline
  --------------------------------------------------------------------------------------------------------------
  b_{i}=v_{i}^{n} + (dt/m_{i})*f_{i} + dt*g_{i}*sum_p(w_{ip})

  Set up vectorXf to contain all b values - Set up one for each direction x, y, z

  Get external force from SimController

  Loop over cell faces - Could use threads here and then split again below

  Loop over all grid cells - Parallel here

    Multiply external force g by e to get in direction of cell face

    Set up sum_p(w)
    Loop over particles in list
      Add to sum_p(w)

    Multiply each part

    Add and set to b for that cell face

  --------------------------------------------------------------------------------------------------------------
  */

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::setUpA_DeviatoricVelocity()
{
  /* Outline
  --------------------------------------------------------------------------------------------------------------
  Loop over grid cells
    Loop over particles list for cell
      Add weight_diff*FE_{p} to dFE_{p} - NB need to set this to zero first in this case

  Loop over grid cells - Not sure if can solve for x, y, z faces separately





  --------------------------------------------------------------------------------------------------------------
  */
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calculate_dR()
{
  /* Outline - For one particle
  ------------------------------------------------------------------------------------------------------------
  Calculate R^T*dF and dF^T*R

  Create b vector

  Create A vector from S from polar decomposition

  Set matrix for R^T*dR

  Solve linear system to obtain x, y and z

  Check these against other three functions, or create second matrix and solve for this too and compare

  Put values into R^T*dR

  Multiply by R to obtain dR

  ------------------------------------------------------------------------------------------------------------
  */
}

//----------------------------------------------------------------------------------------------------------------------
