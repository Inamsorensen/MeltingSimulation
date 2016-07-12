#include "Grid.h"

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcDeviatoricForce(Particle* _particle, Eigen::Vector3f _eVector, Eigen::Vector3f _weightDiff)
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

  //Set final result
  float finalResult;

  //Get deformation gradient elastic, FE, and transpose it
  Eigen::Matrix3f deformationElastic=_particle->getDeformationElastic();
  Eigen::Matrix3f deformationElastic_trans=deformationElastic.transpose();

  //Get particle volume
  float particleVolume=_particle->getVolume();

  //Get deformation gradient deviatoric, J^{-1/d}FE
  Eigen::Matrix3f deformationElastic_Deviatoric=_particle->getDeformationElastic_Deviatoric();

  //Get R from polar decomposition of deformation gradient deviatoric, J^{-1/d}FE
  Eigen::Matrix3f R_deformationElastic_Deviatoric=_particle->getR_deformationElastic_Deviatoric();

  //Get Lame Mu coefficient
  float lameMu=_particle->getLameMu();

  //Calculate dYdFE=2*mu*(FE-RE)
  Eigen::Matrix3f dYdFE=(2.0*lameMu)*(deformationElastic_Deviatoric-R_deformationElastic_Deviatoric);

  //Pass in dYdFE to particle's getZ_DeformEDevDiff. Return is dY^{hat}dFE
  Eigen::Matrix3f dY_hat_dFE=_particle->getZ_DeformEDevDiff(dYdFE);

  //Multiply FE^{T} with differentiated weight
  Eigen::Vector3f result_vec1=deformationElastic_trans*_weightDiff;

  //Multiply result with dY^{hat}dFE
  Eigen::Vector3f result_vec2=dY_hat_dFE*result_vec1;

  //Multiply result with e_{a(i)}
  finalResult=_eVector.dot(result_vec2);

  //Multiply result with particle volume
  finalResult*=particleVolume;

  //Make negative
  finalResult*=(-1.0);

  //Return result
  return finalResult;
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

  //Set e_{a(i)} vectors
  Eigen::Vector3f e_x(1.0, 0.0, 0.0);
  Eigen::Vector3f e_y(0.0, 1.0, 0.0);
  Eigen::Vector3f e_z(0.0, 0.0, 1.0);

  for (int cellIndex=0; cellIndex<(pow(m_noCells, 3)); cellIndex++)
  {
    //Calc total number of cells
    int totNoCells=pow(m_noCells,3);

    //Set up A and b storage
    Eigen::MatrixXf A_X(totNoCells, totNoCells);
    Eigen::VectorXf b_X(totNoCells);
    Eigen::MatrixXf A_Y(totNoCells, totNoCells);
    Eigen::VectorXf b_Y(totNoCells);
    Eigen::MatrixXf A_Z(totNoCells, totNoCells);
    Eigen::VectorXf b_Z(totNoCells);

    //Loop over particles in each face of grid cell
    //Calculate force addition and weight addition

    //Face X

    //Set storage variables
    float forceSumX=0.0;
    float weightSumX=0.0;

    int noParticles_FaceX=m_cellFacesX[cellIndex]->m_interpolationData.size();
    for (int particleIterator=0; particleIterator<noParticles_FaceX; particleIterator++)
    {
      //Get cubic B spline
      float weight=m_cellFacesX[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

      //Get cubic B spline differentiated
      Eigen::Vector3f weight_diff=m_cellFacesX[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline_Diff;

      //Get particle pointer
      Particle* particle=m_cellFacesX[cellIndex]->m_interpolationData[particleIterator]->m_particle;

      //Calculate deviatoric force
      float forceFromParticle=calcDeviatoricForce(particle, e_x, weight_diff);

      //Add force to sum
      forceSumX+=forceFromParticle;

      //Add interpolation weight
      weightSumX+=weight;
    }

    //Face Y
    //Set storage variables
    float forceSumY=0.0;
    float weightSumY=0.0;

    int noParticles_FaceY=m_cellFacesY[cellIndex]->m_interpolationData.size();
    for (int particleIterator=0; particleIterator<noParticles_FaceY; particleIterator++)
    {
      //Get cubic B spline
      float weight=m_cellFacesY[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

      //Get cubic B spline differentiated
      Eigen::Vector3f weight_diff=m_cellFacesY[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline_Diff;

      //Get particle pointer
      Particle* particle=m_cellFacesY[cellIndex]->m_interpolationData[particleIterator]->m_particle;

      //Calculate deviatoric force
      float forceFromParticle=calcDeviatoricForce(particle, e_y, weight_diff);

      //Add force to sum
      forceSumY+=forceFromParticle;

      //Add interpolation weight
      weightSumY+=weight;
    }

    //Face Z
    //Set storage variables
    float forceSumZ=0.0;
    float weightSumZ=0.0;

    int noParticles_FaceZ=m_cellFacesZ[cellIndex]->m_interpolationData.size();
    for (int particleIterator=0; particleIterator<noParticles_FaceZ; particleIterator++)
    {
      //Get cubic B spline
      float weight=m_cellFacesZ[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

      //Get cubic B spline differentiated
      Eigen::Vector3f weight_diff=m_cellFacesZ[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline_Diff;

      //Get particle pointer
      Particle* particle=m_cellFacesZ[cellIndex]->m_interpolationData[particleIterator]->m_particle;

      //Calculate deviatoric force
      float forceFromParticle=calcDeviatoricForce(particle, e_z, weight_diff);

      //Add force to sum
      forceSumZ+=forceFromParticle;

      //Add interpolation weight
      weightSumZ+=weight;
    }


    //Calculate b components for each of the faces
    float velocityX=m_cellFacesX[cellIndex]->m_velocity;
    float massX=m_cellFacesX[cellIndex]->m_mass;
    float velocityY=m_cellFacesY[cellIndex]->m_velocity;
    float massY=m_cellFacesY[cellIndex]->m_mass;
    float velocityZ=m_cellFacesZ[cellIndex]->m_velocity;
    float massZ=m_cellFacesZ[cellIndex]->m_mass;

    float bComponentX=calcBComponent_DeviatoricVelocity(velocityX, massX, forceSumX, weightSumX, e_x);
    float bComponentY=calcBComponent_DeviatoricVelocity(velocityY, massY, forceSumY, weightSumY, e_y);
    float bComponentZ=calcBComponent_DeviatoricVelocity(velocityZ, massZ, forceSumZ, weightSumZ, e_z);

    //Insert into b vector
    b_X(cellIndex)=bComponentX;
    b_Y(cellIndex)=bComponentY;
    b_Z(cellIndex)=bComponentZ;

  }


}

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcBComponent_DeviatoricVelocity(float _velocity, float _mass, float _deviatoricForce, float _sumWeight, Eigen::Vector3f _eVector)
{
  /* Outline
  --------------------------------------------------------------------------------------------------------------
  b_{i}=v_{i}^{n} + (dt/m_{i})*f_{i} + dt*g_{i}*sum_p(w_{ip})
  --------------------------------------------------------------------------------------------------------------
  */

  //Set result to return
  float result;

  //External forces component
  float externalForce=m_externalForce.dot(_eVector);
  externalForce*=(m_dt*_sumWeight);

  //Deviatoric force component
  float deviatoricForce=(m_dt/_mass)*_deviatoricForce;

  //Calculate b component
  result=_velocity+deviatoricForce+externalForce;

  //Return result
  return result;
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
