#include "Grid.h"

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcDeviatoricForce(Particle* _particle)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  f_{i}=-sum_p(V_{p}*e_{a(i)}^T*dYdFE_{p}:B_{kmij}*FE_{p}^T*cubicBSpline_Diff_{ip}

  dYdFE_{p}=2*lame_mu_{p}*JE^a*FE_{p} - 2*lame_mu_{p}*RE_{p} where a=-1/d where d=3? for dimensions and RE is the polar decomposition of JE^a*FE

  B_{kmij}=JE^a*I + a*JE^a*FE_{p}^-T*FE_{p} - Stored in particle

  Loop over all grid cells
    Loop over all particles in grid cell

      calc dYdFE_{p}

      special multiply dYdFE_{p} and deformationElastic_Deviatoric_Diff

      multiply with FE^T

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
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Need some of the same variables for calculating deviatoric force and elements of A for implicit velocity calc
  Hence this function loops over all particles inside each grid cell and all grid cells. Other functions accessed
  on a per particle basis

  Loop over grid cells

    Need to set up A and b storage

    Loop over particles in grid cell

      Run calc deviatoric force

      Calc sum weights for external forces

    End loop particles

    Run calc B element

    Input into B vector - one for each face

  End loop grid cells

  Set up A matrix if implicit calc on

    Loop over grid cells
      Loop over grid cells
        Run calc A element
        Input into A vector

  Run implicit solve - one thread for each face?

  Store results
    Loop over grid cells again to store

  ----------------------------------------------------------------------------------------------------------------
  */


}

//----------------------------------------------------------------------------------------------------------------------

void Grid::setUpB_DeviatoricVelocity()
{
  /* Outline
  --------------------------------------------------------------------------------------------------------------
  b_{i}=v_{i}^{n} + (dt/m_{i})*f_{i} + dt*g_{i}*sum_p(w_{ip})

  Get external force from SimController

  Loop over cell faces - Could use threads here and then split again below

    Multiply external force g by e to get in direction of cell face and my sum weights

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
  For each face

    calculate dF=B:Z

    calculate_dR();

    calculate d2YdFE2

    calculate C_hat:Z by calculating each part and adding

    Multiply by FE^T and weight_diff

    Multiply by e_a^T and V_{p}

    Add to A element for face

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
