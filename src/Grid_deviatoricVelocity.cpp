#include "Grid.h"


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


  TODO: Set boundary velocities to go into b

  ----------------------------------------------------------------------------------------------------------------
  */

  //Set e_{a(i)} vectors
  Eigen::Vector3f e_x(1.0, 0.0, 0.0);
  Eigen::Vector3f e_y(0.0, 1.0, 0.0);
  Eigen::Vector3f e_z(0.0, 0.0, 1.0);

  //Set up A and b storage
  Eigen::VectorXf b_X(m_totNoCells);
  Eigen::VectorXf b_Y(m_totNoCells);
  Eigen::VectorXf b_Z(m_totNoCells);

  Eigen::MatrixXf A_X(m_totNoCells, m_totNoCells);
  Eigen::MatrixXf A_Y(m_totNoCells, m_totNoCells);
  Eigen::MatrixXf A_Z(m_totNoCells, m_totNoCells);

  //Initialise all to zero
  b_X.setZero();
  b_Y.setZero();
  b_Z.setZero();
  A_X.setZero();
  A_Y.setZero();
  A_Z.setZero();

  //Set whether implicit or explicit integration
  bool implicitUpdate=false;
//  bool implicitUpdate=true;

//#pragma omp parallel for
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    //Test parallel
//    printf("The parallel region is executed by thread %d\n", omp_get_thread_num());
//    std::cout<<"test\n";


    //Store velocities as step n before update them
    m_cellFacesX[cellIndex]->m_previousVelocity=m_cellFacesX[cellIndex]->m_velocity;
    m_cellFacesY[cellIndex]->m_previousVelocity=m_cellFacesY[cellIndex]->m_velocity;
    m_cellFacesZ[cellIndex]->m_previousVelocity=m_cellFacesZ[cellIndex]->m_velocity;

//    //Calculate number of particles in faces of cellIndex
//    int noParticles_FaceX=m_cellFacesX[cellIndex]->m_interpolationData.size();
//    int noParticles_FaceY=m_cellFacesY[cellIndex]->m_interpolationData.size();
//    int noParticles_FaceZ=m_cellFacesZ[cellIndex]->m_interpolationData.size();

    //Get state of cell face instead
    State state_FaceX=m_cellFacesX[cellIndex]->m_state;
    State state_FaceY=m_cellFacesY[cellIndex]->m_state;
    State state_FaceZ=m_cellFacesZ[cellIndex]->m_state;

    //Only set B components for non-empty cells
//    if (noParticles_FaceX!=0 && noParticles_FaceY!=0 && noParticles_FaceZ!=0)
    if (state_FaceX==State::Interior)
    {
      b_X(cellIndex)=calcBComponent_DeviatoricVelocity(m_cellFacesX[cellIndex], e_x);
    }

    if (state_FaceY==State::Interior)
    {
      b_Y(cellIndex)=calcBComponent_DeviatoricVelocity(m_cellFacesY[cellIndex], e_y);
    }

    if (state_FaceZ==State::Interior)
    {
      b_Z(cellIndex)=calcBComponent_DeviatoricVelocity(m_cellFacesZ[cellIndex], e_z);
    }

//    if (state_FaceX!=State::Empty || state_FaceY!=State::Empty || state_FaceZ!=State::Empty)
//    {
//      B_components=calcBComponent_DeviatoricVelocity(cellIndex);

//      //Insert into b vector
//      b_X(cellIndex)=B_components(0);
//      b_Y(cellIndex)=B_components(1);
//      b_Z(cellIndex)=B_components(2);
//    }


    //Get A components for implicit update
    if (implicitUpdate==true)
    {
      //Calculate number of particles in faces of cellIndex
      int noParticles_FaceX=m_cellFacesX[cellIndex]->m_interpolationData.size();
      int noParticles_FaceY=m_cellFacesY[cellIndex]->m_interpolationData.size();
      int noParticles_FaceZ=m_cellFacesZ[cellIndex]->m_interpolationData.size();

      //Only set A components for non-empty cells
//      if (noParticles_FaceX!=0 && noParticles_FaceY!=0 && noParticles_FaceZ!=0)
      if (state_FaceX!=State::Empty || state_FaceY!=State::Empty || state_FaceZ!=State::Empty)
      {
        calcAComponent_DeviatoricVelocity(cellIndex, noParticles_FaceX, noParticles_FaceY, noParticles_FaceZ, A_X, A_Y, A_Z);
      }
    }

    else
    {
//      float massX=m_cellFacesX[cellIndex]->m_mass;
//      float massY=m_cellFacesY[cellIndex]->m_mass;
//      float massZ=m_cellFacesZ[cellIndex]->m_mass;

      float velocityX=0.0;
      float velocityY=0.0;
      float velocityZ=0.0;

//      if (massX!=0)
      if (state_FaceX==State::Interior)
      {
        velocityX=(b_X(cellIndex)/m_cellFacesX[cellIndex]->m_mass);
      }

//      if (massY!=0)
      if (state_FaceY==State::Interior)
      {
        velocityY=(b_Y(cellIndex)/m_cellFacesY[cellIndex]->m_mass);
      }

//      if (massZ!=0)
      if (state_FaceZ==State::Interior)
      {
        velocityZ=(b_Z(cellIndex)/m_cellFacesZ[cellIndex]->m_mass);
      }

      explicitUpdateVelocity(cellIndex, velocityX, velocityY, velocityZ);

    }

  }

  if (implicitUpdate==true)
  {
    implicitUpdateVelocity(A_X, b_X, A_Y, b_Y, A_Z, b_Z);

  }


}



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

  //Multiply FE^{T} with differentiated weight
  Eigen::Vector3f result_vec1=deformationElastic_trans*_weightDiff;

  //Get dY^{hat}dFE from particle
  Eigen::Matrix3f dY_hat_dFE=_particle->getPotentialEnergyDiff();

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


//----------------------------------------------------------------------------------------------------------------------

float Grid::calcBComponent_DeviatoricVelocity(CellFace *_cellFace, Eigen::Vector3f _eVector)
{
  /* Outline
  --------------------------------------------------------------------------------------------------------------
  b_{i}=v_{i}^{n} + (dt/m_{i})*f_{i} + dt*g_{i}*sum_p(w_{ip})

  where f_{i} is the force in the i face of direction X/Y/Z. Calculated by calcDeviatoricForce

  Returns vector where B_components(0) is the B component for the X direction etc.
  --------------------------------------------------------------------------------------------------------------
  */

  //Components to return
  float BComponent=0.0;


  //Set storage variables
  float forceSum=0.0;
  float weightSum=0.0;

  int noParticles=_cellFace->m_interpolationData.size();
  for (int particleIterator=0; particleIterator<noParticles; particleIterator++)
  {
    //Get cubic B spline
    float weight=_cellFace->m_interpolationData[particleIterator]->m_cubicBSpline;

    //Get cubic B spline differentiated
    Eigen::Vector3f weight_diff=_cellFace->m_interpolationData[particleIterator]->m_cubicBSpline_Diff;

    //Get particle pointer
    Particle* particle=_cellFace->m_interpolationData[particleIterator]->m_particle;

    float forceFromParticle=calcDeviatoricForce(particle, _eVector, weight_diff);

    //Add force to sum
    forceSum+=forceFromParticle;

    //Add interpolation weight
    weightSum+=weight;
  }

  //Determine if cell is empty, if so B_component is zero. Need to do this in case collision cell is empty
  ///Might be able to take this out since should be no particles inside collision cells
  if (noParticles!=0)
  {
    float cellMass=_cellFace->m_mass;

    //External forces component
    float externalForce=0.0;
    externalForce=m_externalForce.dot(_eVector);
//    externalForceX*=(m_dt*weightSumX);
    externalForce*=(m_dt*weightSum*cellMass);

    //Deviatoric force component
//    float deviatoricForceX=(m_dt/cellMassX)*forceSumX;
    float deviatoricForce=(m_dt*forceSum);

    //Calculate b component
//    B_components(0)=m_cellFacesX[_cellIndex]->m_velocity+deviatoricForceX+externalForceX;
    BComponent=(_cellFace->m_velocity*cellMass) + deviatoricForce + externalForce;
  }


//  //Face Y
//  //Set storage variables
//  float forceSumY=0.0;
//  float weightSumY=0.0;

//  int noParticles_FaceY=m_cellFacesY[_cellIndex]->m_interpolationData.size();
//  for (int particleIterator=0; particleIterator<noParticles_FaceY; particleIterator++)
//  {
//    //Get cubic B spline
//    float weight=m_cellFacesY[_cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

//    //Get cubic B spline differentiated
//    Eigen::Vector3f weight_diff=m_cellFacesY[_cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline_Diff;

//    //Get particle pointer
//    Particle* particle=m_cellFacesY[_cellIndex]->m_interpolationData[particleIterator]->m_particle;

//    float forceFromParticle=calcDeviatoricForce(particle, e_y, weight_diff);

//    //Add force to sum
//    forceSumY+=forceFromParticle;

//    //Add interpolation weight
//    weightSumY+=weight;
//  }

//  //Determine if cell is empty, if so B_component is zero
//  if (noParticles_FaceY!=0)
//  {
//    float cellMassY=m_cellFacesY[_cellIndex]->m_mass;

//    //External forces component
//    float externalForceY=0.0;
//    externalForceY=m_externalForce.dot(e_y);
////    externalForceY*=(m_dt*weightSumY);
//    externalForceY*=(m_dt*weightSumY*cellMassY);

//    //Deviatoric force component
////    float deviatoricForceY=(m_dt/cellMassY)*forceSumY;
//    float deviatoricForceY=(m_dt*forceSumY);

//    //Calculate b component
////    B_components(1)=m_cellFacesY[_cellIndex]->m_velocity + deviatoricForceY + externalForceY;
//    B_components(1)=(m_cellFacesY[_cellIndex]->m_velocity*cellMassY) + deviatoricForceY + externalForceY;
//  }



//  //Face Z
//  //Set storage variables
//  float forceSumZ=0.0;
//  float weightSumZ=0.0;

//  int noParticles_FaceZ=m_cellFacesZ[_cellIndex]->m_interpolationData.size();
//  for (int particleIterator=0; particleIterator<noParticles_FaceZ; particleIterator++)
//  {
//    //Get cubic B spline
//    float weight=m_cellFacesZ[_cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

//    //Get cubic B spline differentiated
//    Eigen::Vector3f weight_diff=m_cellFacesZ[_cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline_Diff;

//    //Get particle pointer
//    Particle* particle=m_cellFacesZ[_cellIndex]->m_interpolationData[particleIterator]->m_particle;

//    float forceFromParticle=calcDeviatoricForce(particle, e_z, weight_diff);

//    //Add force to sum
//    forceSumZ+=forceFromParticle;

//    //Add interpolation weight
//    weightSumZ+=weight;
//  }

//  //Determine if cell is empty, if so B_component is zero
//  if (noParticles_FaceZ!=0)
//  {
//    float cellMassZ=m_cellFacesZ[_cellIndex]->m_mass;

//    //External forces component
//    float externalForceZ=0.0;
//    externalForceZ=m_externalForce.dot(e_z);
////    externalForceZ*=(m_dt*weightSumZ);
//    externalForceZ*=(m_dt*weightSumZ*cellMassZ);

//    //Deviatoric force component
////    float deviatoricForceZ=(m_dt/cellMassZ)*forceSumZ;
//    float deviatoricForceZ=(m_dt*forceSumZ);

//    //Calculate b component
////    B_components(2)=m_cellFacesZ[_cellIndex]->m_velocity + deviatoricForceZ + externalForceZ;
//    B_components(2)=(m_cellFacesZ[_cellIndex]->m_velocity*cellMassZ) + deviatoricForceZ + externalForceZ;
//  }


  //Return result
  return BComponent;
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcAComponent_DeviatoricVelocity(int _cellIndex, int _noParticlesFaceX, int _noParticlesFaceY, int _noParticlesFaceZ, Eigen::MatrixXf &o_AX, Eigen::MatrixXf &o_AY, Eigen::MatrixXf &o_AZ)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------------
  Loops through the neighbours at +-2 increments from given cellIndex in all directions

  Check that cell neighbour isn't empty

  If not empty, loops through particles for both faces to find the ones that are the same
    NB! Guessing this is slow. Is there a faster way of doing this?

  Once same particle in both is found, calculate contribution from that particle

  ------------------------------------------------------------------------------------------------------------
  */

  //Set e_{a(i)} vectors
  Eigen::Vector3f e_x(1.0, 0.0, 0.0);
  Eigen::Vector3f e_y(0.0, 1.0, 0.0);
  Eigen::Vector3f e_z(0.0, 0.0, 1.0);

  //Get indices of current cell
  int iIndex=m_cellCentres[_cellIndex]->m_iIndex;
  int jIndex=m_cellCentres[_cellIndex]->m_jIndex;
  int kIndex=m_cellCentres[_cellIndex]->m_kIndex;

  //Set min and max increments for ijk
  int iMinIncrement=-2;
  int jMinIncrement=-2;
  int kMinIncrement=-2;
  int iMaxIncrement=2;
  int jMaxIncrement=2;
  int kMaxIncrement=2;

  //Make sure doesn't loop over outside cells.
  //Sort of makes the non-existing cells seem like collision cells as will add nothing to volume
  // i direction
  if (iIndex==0)
  {
    iMinIncrement=0;
  }
  if (iIndex==1)
  {
    iMinIncrement=-1;
  }
  if (iIndex==(m_noCells-2))
  {
    iMaxIncrement=1;
  }
  if (iIndex==(m_noCells-1))
  {
    iMaxIncrement=0;
  }

  //j direction
  if (jIndex==0)
  {
    jMinIncrement=0;
  }
  if (jIndex==1)
  {
    jMinIncrement=-1;
  }
  if (jIndex==(m_noCells-2))
  {
    jMaxIncrement=1;
  }
  if (jIndex==(m_noCells-1))
  {
    jMaxIncrement=0;
  }

  // k direction
  if (kIndex==0)
  {
    kMinIncrement=0;
  }
  if (kIndex==1)
  {
    kMinIncrement=-1;
  }
  if (kIndex==(m_noCells-2))
  {
    kMaxIncrement=1;
  }
  if (kIndex==(m_noCells-1))
  {
    kMaxIncrement=0;
  }

  //Get state of cell faces
  State state_FaceX=m_cellFacesX[_cellIndex]->m_state;
  State state_FaceY=m_cellFacesY[_cellIndex]->m_state;
  State state_FaceZ=m_cellFacesZ[_cellIndex]->m_state;


  //Loop over neighbouring cells that can have same particles in them
  for (int kIndexIncrement=kMinIncrement; kIndexIncrement<(kMaxIncrement+1); kIndexIncrement++)
  {
    for (int jIndexIncrement=jMinIncrement; jIndexIncrement<(jMaxIncrement+1); jIndexIncrement++)
    {
      for (int iIndexIncrement=iMinIncrement; iIndexIncrement<(iMaxIncrement+1); iIndexIncrement++)
      {
        //Initialise A components
        float AcomponentX=0.0;
        float AcomponentY=0.0;
        float AcomponentZ=0.0;

        //Get index of neighbour
        int neighbourCellIndex=MathFunctions::getVectorIndex(iIndex+iIndexIncrement, jIndex+jIndexIncrement, kIndex+kIndexIncrement, m_noCells);

        //Get state of neighbour cell
        State state_FaceX_neighbour=m_cellFacesX[neighbourCellIndex]->m_state;
        State state_FaceY_neighbour=m_cellFacesY[neighbourCellIndex]->m_state;
        State state_FaceZ_neighbour=m_cellFacesZ[neighbourCellIndex]->m_state;

        //Only insert an A component if neighbour face isn't colliding
        if (state_FaceX==State::Interior && state_FaceX_neighbour==State::Interior)
        {
//          int noParticlesFaceX_neighbour=m_cellFacesX[neighbourCellIndex]->m_interpolationData.size();

          //Face X
//          if (noParticlesFaceX_neighbour!=0 && _noParticlesFaceX!=0)
//          if (state_FaceX==State::Interior && state_FaceX_neighbour==State::Interior)
//          {

            //Find same particle in list
            for (int particleIterator_i=0; particleIterator_i<_noParticlesFaceX; particleIterator_i++)
            {
              //Get id of particle i
              unsigned int particleId_i=m_cellFacesX[_cellIndex]->m_interpolationData[particleIterator_i]->m_particle->getId();

              bool isFound=false;
              unsigned int particleId_j;
              searchCellsForCommonParticle(particleId_i, m_cellFacesX[neighbourCellIndex], particleId_j, isFound);

              if (isFound==true)
              {
                //Get weights and mass
                Eigen::Vector3f weight_i_diff=m_cellFacesX[_cellIndex]->m_interpolationData[particleIterator_i]->m_cubicBSpline_Diff;
                Eigen::Vector3f weight_j_diff=m_cellFacesX[neighbourCellIndex]->m_interpolationData[particleId_j]->m_cubicBSpline_Diff;

                //Get particle pointer
                Particle* commonParticle=m_cellFacesX[_cellIndex]->m_interpolationData[particleIterator_i]->m_particle;

                //Calculate A value for specific particle
                float AValue_particle=calcAValue_DeviatoricVelocity(commonParticle, weight_i_diff, weight_j_diff, e_x);

                //Add to A_ij value
                AcomponentX+=AValue_particle;
              }
            }

            //Add mass to diagonal elements
            if (_cellIndex==neighbourCellIndex)
            {
              float mass_i=m_cellFacesX[_cellIndex]->m_mass;
              AcomponentX+=mass_i;
            }

            //Insert A component to matrix
            o_AX(_cellIndex, neighbourCellIndex)=AcomponentX;
//          }
        }

        //Face Y
        if (state_FaceY==State::Interior && state_FaceY_neighbour==State::Interior)
        {
//          int noParticlesFaceY_neighbour=m_cellFacesY[neighbourCellIndex]->m_interpolationData.size();

//          if (noParticlesFaceY_neighbour!=0 && _noParticlesFaceY!=0)
//          {
            //Find same particle in list
            for (int particleIterator_i=0; particleIterator_i<_noParticlesFaceY; particleIterator_i++)
            {
              //Get id of particle i
              unsigned int particleId_i=m_cellFacesY[_cellIndex]->m_interpolationData[particleIterator_i]->m_particle->getId();

              bool isFound=false;
              unsigned int particleId_j;
              searchCellsForCommonParticle(particleId_i, m_cellFacesY[neighbourCellIndex], particleId_j, isFound);

              if (isFound==true)
              {
                //Get weights and mass
                Eigen::Vector3f weight_i_diff=m_cellFacesY[_cellIndex]->m_interpolationData[particleIterator_i]->m_cubicBSpline_Diff;
                Eigen::Vector3f weight_j_diff=m_cellFacesY[neighbourCellIndex]->m_interpolationData[particleId_j]->m_cubicBSpline_Diff;

                //Get particle pointer
                Particle* commonParticle=m_cellFacesY[_cellIndex]->m_interpolationData[particleIterator_i]->m_particle;

                //Calculate A value for specific particle
                float AValue_particle=calcAValue_DeviatoricVelocity(commonParticle, weight_i_diff, weight_j_diff, e_y);

                //Add to A_ij value
                AcomponentY+=AValue_particle;
              }
            }

            //Add mass to diagonal elements
            if (_cellIndex==neighbourCellIndex)
            {
              float mass_i=m_cellFacesY[_cellIndex]->m_mass;
              AcomponentY+=mass_i;
            }

            //Insert A component to matrix
            o_AY(_cellIndex, neighbourCellIndex)=AcomponentY;
//          }

        }


        //Face Z
        if (state_FaceZ==State::Interior && state_FaceZ_neighbour==State::Interior)
        {
//          int noParticlesFaceZ_neighbour=m_cellFacesZ[neighbourCellIndex]->m_interpolationData.size();

//          if (noParticlesFaceZ_neighbour!=0 && _noParticlesFaceZ!=0)
//          {
            //Find same particle in list
            for (int particleIterator_i=0; particleIterator_i<_noParticlesFaceZ; particleIterator_i++)
            {
              //Get id of particle i
              unsigned int particleId_i=m_cellFacesZ[_cellIndex]->m_interpolationData[particleIterator_i]->m_particle->getId();

              bool isFound=false;
              unsigned int particleId_j;
              searchCellsForCommonParticle(particleId_i, m_cellFacesZ[neighbourCellIndex], particleId_j, isFound);

              if (isFound==true)
              {
                //Get weights and mass
                Eigen::Vector3f weight_i_diff=m_cellFacesZ[_cellIndex]->m_interpolationData[particleIterator_i]->m_cubicBSpline_Diff;
                Eigen::Vector3f weight_j_diff=m_cellFacesZ[neighbourCellIndex]->m_interpolationData[particleId_j]->m_cubicBSpline_Diff;

                //Get particle pointer
                Particle* commonParticle=m_cellFacesZ[_cellIndex]->m_interpolationData[particleIterator_i]->m_particle;

                //Calculate A value for specific particle
                float AValue_particle=calcAValue_DeviatoricVelocity(commonParticle, weight_i_diff, weight_j_diff, e_z);

                //Add to A_ij value
                AcomponentZ+=AValue_particle;
              }
            }

            //Add mass to diagonal elements
            if (_cellIndex==neighbourCellIndex)
            {
              float mass_i=m_cellFacesZ[_cellIndex]->m_mass;
              AcomponentZ+=mass_i;
            }

            //Insert A component to matrix
            o_AZ(_cellIndex, neighbourCellIndex)=AcomponentZ;
//          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcAValue_DeviatoricVelocity(Particle *_particle, Eigen::Vector3f _weight_i_diff, Eigen::Vector3f _weight_j_diff, Eigen::Vector3f _eVector)
{
  /* Outline
  --------------------------------------------------------------------------------------------------------------
  For each face

    calculate Z

    calculate deltaJF=B:Z

    calculate_dR();

    calculate d2YdFE2

    calculate C_hat:Z by calculating each part and adding

    Multiply by FE^T and weight_diff

    Multiply by e_a^T and V_{p}

    Add to A element for face

  --------------------------------------------------------------------------------------------------------------
  */

  //Set up storage for different parts
  Eigen::Matrix3f part1;
  Eigen::Matrix3f part2;
  Eigen::Matrix3f part3;
  Eigen::Matrix3f part4;

  //Get parameters from particle
  Eigen::Matrix3f deformationGradElastic=_particle->getDeformationElastic();
  Eigen::Matrix3f R_deformElastic_Deviatoric=_particle->getR_deformationElastic_Deviatoric();
  Eigen::Matrix3f S_deformElastic_Deviatoric=_particle->getS_deformationElastic_Deviatoric();
  float lameMu=_particle->getLameMu();
  float dimension=_particle->getDimension();
  float detDeformElastic=_particle->getDetDeformationElastic();

  Eigen::Matrix3f deformElastic_trans=deformationGradElastic.transpose();
  Eigen::Matrix3f deformElastic_trans_inverse=deformElastic_trans.inverse();

  //Calculate Z
  Eigen::Matrix3f Z_matrix;
  Z_matrix(0,0)=_eVector(0)*_weight_j_diff(0);
  Z_matrix(0,1)=_eVector(0)*_weight_j_diff(1);
  Z_matrix(0,2)=_eVector(0)*_weight_j_diff(2);
  Z_matrix(1,0)=_eVector(1)*_weight_j_diff(0);
  Z_matrix(1,1)=_eVector(1)*_weight_j_diff(1);
  Z_matrix(1,2)=_eVector(1)*_weight_j_diff(2);
  Z_matrix(2,0)=_eVector(2)*_weight_j_diff(0);
  Z_matrix(2,1)=_eVector(2)*_weight_j_diff(1);
  Z_matrix(2,2)=_eVector(2)*_weight_j_diff(2);

  Z_matrix=(Z_matrix*deformationGradElastic);

  Eigen::Matrix3f Z_matrix_trans=Z_matrix.transpose();


  //Calculate deltaJF=B:Z_matrix
  Eigen::Matrix3f deltaJF=_particle->getDeformEDevDiff_Z(Z_matrix);

  //Calculate first part of differential
  //Calculate deltaR
  Eigen::Matrix3f deltaR=calculate_dR(deltaJF, R_deformElastic_Deviatoric, S_deformElastic_Deviatoric);

  //Calculate d2YdF2
  Eigen::Matrix3f d2YdF2=(2.0*lameMu)*(deltaJF-deltaR);

  //Find part 1 of differential
  part1=_particle->getZ_DeformEDevDiff(d2YdF2);


  //Get dYdF which is used in part 2,3,4
//  Eigen::Matrix3f dYdFE=(2.0*lameMu)*(deformationElastic_Deviatoric-R_deformElastic_Deviatoric);
  Eigen::Matrix3f dYdFE=_particle->getPotentialEnergyDiff();

  //Calculate a=-1/dimensions and J^a
  float aConstant=(-1.0/dimension);
  float JaConstant=pow(detDeformElastic, aConstant);

  //Calculate part 2
  Eigen::Matrix3f part2_1=MathFunctions::matrixElementMultiplication(deformElastic_trans_inverse, Z_matrix);
  Eigen::Matrix3f part2_2=_particle->getZ_DeformEDevDiff(dYdFE);
  part2=part2_1*part2_2;
  part2*=aConstant;

  //Calculate part 3
  part3=MathFunctions::matrixElementMultiplication(dYdFE, Z_matrix);
  part3=part3*deformElastic_trans_inverse;
  part3*=(aConstant*JaConstant);

  //Calculate part 4
  part4=MathFunctions::matrixElementMultiplication(dYdFE, deformationGradElastic);
  part4=part4*deformElastic_trans_inverse;
  part4=part4*Z_matrix_trans;
  part4=part4*deformElastic_trans_inverse;
  part4*=(aConstant*JaConstant);
  part4*=(-1.0);


  //Add all parts to find C_hat:Z
  Eigen::Matrix3f Ap_matrix=(part1+part2+part3+part4);


  //Multiply Ap matrix with other particle dependent variables
  Eigen::Matrix3f Ap_deformElastTrans=Ap_matrix*deformElastic_trans;

  Eigen::Vector3f Ap_deformElastTrans_weightDiff=Ap_deformElastTrans*_weight_i_diff;

  float Acomponent=_eVector.dot(Ap_deformElastTrans_weightDiff);

  Acomponent*=_particle->getVolume();

//  //Multiply by deltaT^2/2*mass_i
//  Acomponent*=((pow(m_dt,2.0))/(2.0*_mass_i));
  //Multiply by deltaT^2/2
  Acomponent*=((pow(m_dt,2.0))/2.0);

  //Return result
  return Acomponent;


}

//----------------------------------------------------------------------------------------------------------------------

Eigen::Matrix3f Grid::calculate_dR(const Eigen::Matrix3f &_deltaDeformElastic_Deviatoric, const Eigen::Matrix3f &_R_deformElastic_Deviatoric, const Eigen::Matrix3f &_S_deformElastic_Deviatoric)
{
  /* Outline - For one particle
  ------------------------------------------------------------------------------------------------------------
  Calculate R^T*dJF and dJF^T*R

  Create b vector

  Create A vector from S from polar decomposition

  Set matrix for R^T*dR

  Verify that R^T*dJF-dJF^T*R is anti symmetric. Otherwise the remaining 3 equations are not linearly dependent

  Solve linear system to obtain x, y and z

  Check these against other three functions, or create second matrix and solve for this too and compare

  Put values into R^T*dR

  Multiply by R to obtain dR

  ------------------------------------------------------------------------------------------------------------
  */

  //Calculate R^T*dJF and dJF^T*R
  Eigen::Matrix3f R_transpose=_R_deformElastic_Deviatoric.transpose();
  Eigen::Matrix3f deltaJF_transpose=_deltaDeformElastic_Deviatoric.transpose();

  Eigen::Matrix3f B_matrix=(R_transpose*_deltaDeformElastic_Deviatoric);
  B_matrix-=(deltaJF_transpose*_R_deformElastic_Deviatoric);

  //Create A matrix and b vector
  Eigen::Matrix3f A_matrix;
  Eigen::Vector3f b_vector;
  Eigen::Vector3f solution_vector;

  //Insert elements to A
  ///NB! Double check these elements
  A_matrix(0,0)=_S_deformElastic_Deviatoric(0,0)+_S_deformElastic_Deviatoric(1,1);
  A_matrix(0,1)=_S_deformElastic_Deviatoric(1,2);
  A_matrix(0,2)=(-1.0)*_S_deformElastic_Deviatoric(0,2);
  A_matrix(1,0)=A_matrix(0,1);
  A_matrix(1,1)=_S_deformElastic_Deviatoric(0,0)+_S_deformElastic_Deviatoric(2,2);
  A_matrix(1,2)=_S_deformElastic_Deviatoric(0,1);
  A_matrix(2,0)=A_matrix(0,2);
  A_matrix(2,1)=A_matrix(1,2);
  A_matrix(2,2)=_S_deformElastic_Deviatoric(1,1)+_S_deformElastic_Deviatoric(2,2);

  //Insert elements to b vector
  b_vector(0)=B_matrix(0,1);
  b_vector(1)=B_matrix(0,2);
  b_vector(2)=B_matrix(1,2);

  //Initialise solution to zero
  solution_vector.setZero();

  //Verify that B matrix is antisymmetric. Otherwise non-dependent equations.
  //This check might not work due to floating point inaccuracies

  float tolerance=pow(10.0,-7.0);

  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      //Calculate B(i,j)+B(j,i). Should be zero cause antisymmetric matrices have opposite negatives
      //and zero on diagonal
      float antisymmetryCheck=B_matrix(i,j)+B_matrix(j,i);

      if (antisymmetryCheck>tolerance)
      {
        throw std::invalid_argument("The vector used to calculate delta R is not anti symmetric.");
      }
    }
  }

  //Use linear solver to find x, y, z for RT*deltaR
  MathFunctions::linearSystemSolve(A_matrix, b_vector, solution_vector);


  //Insert solution to RT*deltaR
  Eigen::Matrix3f Rtrans_deltaR;
  Rtrans_deltaR(0,0)=0.0;
  Rtrans_deltaR(0,1)=solution_vector(0);
  Rtrans_deltaR(0,2)=solution_vector(1);
  Rtrans_deltaR(1,0)=(-1.0)*Rtrans_deltaR(0,1);
  Rtrans_deltaR(1,1)=0.0;
  Rtrans_deltaR(1,2)=solution_vector(2);
  Rtrans_deltaR(2,0)=(-1.0)*Rtrans_deltaR(0,2);
  Rtrans_deltaR(2,1)=(-1.0)*Rtrans_deltaR(1,2);
  Rtrans_deltaR(2,2)=0.0;

  //Multiply with R to get deltaR
  Eigen::Matrix3f deltaR=_R_deformElastic_Deviatoric*Rtrans_deltaR;

  //Return result
  return deltaR;
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::explicitUpdateVelocity(int _cellIndex, float _velocityX, float _velocityY, float _velocityZ)
{
  m_cellFacesX[_cellIndex]->m_velocity=_velocityX;
  m_cellFacesY[_cellIndex]->m_velocity=_velocityY;
  m_cellFacesZ[_cellIndex]->m_velocity=_velocityZ;
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::implicitUpdateVelocity(const Eigen::MatrixXf &_A_X, const Eigen::VectorXf &_bVector_X, const Eigen::MatrixXf &_A_Y, const Eigen::VectorXf &_bVector_Y, const Eigen::MatrixXf &_A_Z, const Eigen::VectorXf &_bVector_Z)
{
//  //Set up face normals
//  Eigen::Vector3f e_x(1.0, 0.0, 0.0);
//  Eigen::Vector3f e_y(0.0, 1.0, 0.0);
//  Eigen::Vector3f e_z(0.0, 0.0, 1.0);

//  //Set up A storage
//  Eigen::MatrixXf A_X(m_totNoCells, m_totNoCells);
//  Eigen::MatrixXf A_Y(m_totNoCells, m_totNoCells);
//  Eigen::MatrixXf A_Z(m_totNoCells, m_totNoCells);
//  A_X.setZero();
//  A_Y.setZero();
//  A_Z.setZero();

////  printf("Nested parallelism is %s\n", omp_get_nested() ? "supported" : "not supported");

//  omp_set_nested(1);
//#pragma omp parallel for
//  for (int cellIndex_i=0; cellIndex_i<m_totNoCells; cellIndex_i++)
//  {
//    //Test parallel
////    printf("Thread %d executes outer parallel region\n", omp_get_thread_num());

//    //Calculate number of particles in faces of cellIndex_j
//    int noParticles_FaceX_i=m_cellFacesX[cellIndex_i]->m_interpolationData.size();
//    int noParticles_FaceY_i=m_cellFacesY[cellIndex_i]->m_interpolationData.size();
//    int noParticles_FaceZ_i=m_cellFacesZ[cellIndex_i]->m_interpolationData.size();

//    //If all faces are empty then skip
//    if (noParticles_FaceX_i==0 && noParticles_FaceY_i==0 && noParticles_FaceZ_i==0)
//    {
//      continue;
//    }

//    //Loop over cells again to calculate A components
//#pragma omp parallel for
//    for (int cellIndex_j=0; cellIndex_j<m_totNoCells; cellIndex_j++)
//    {
//      //Test parallel
////      printf("Thread %d executes inner parallel region\n", omp_get_thread_num());

//      //Calculate number of particles in faces of cellIndex_j
//      int noParticles_FaceX_j=m_cellFacesX[cellIndex_j]->m_interpolationData.size();
//      int noParticles_FaceY_j=m_cellFacesY[cellIndex_j]->m_interpolationData.size();
//      int noParticles_FaceZ_j=m_cellFacesZ[cellIndex_j]->m_interpolationData.size();

//      //Face X
//      //Loop over particles in cell face X for cellIndex
//      for (int particleIterator_i=0; particleIterator_i<noParticles_FaceX_i; particleIterator_i++)
//      {
//        for (int particleIterator_j=0; particleIterator_j<noParticles_FaceX_j; particleIterator_j++)
//        {
//          Particle* particleX_i=m_cellFacesX[cellIndex_i]->m_interpolationData[particleIterator_i]->m_particle;
//          Particle* particleX_j=m_cellFacesX[cellIndex_j]->m_interpolationData[particleIterator_j]->m_particle;

//          if (particleX_i==particleX_j)
//          {
//            //Calculate Aij component from p

//            //Get interpolation weights
//            Eigen::Vector3f cubicBSplineDiff_i=m_cellFacesX[cellIndex_i]->m_interpolationData[particleIterator_i]->m_cubicBSpline_Diff;
//            Eigen::Vector3f cubicBSplineDiff_j=m_cellFacesX[cellIndex_j]->m_interpolationData[particleIterator_j]->m_cubicBSpline_Diff;

//            float mass_faceX_i=m_cellFacesX[cellIndex_i]->m_mass;
//            float Acomponent=calcAComponent_DeviatoricVelocity(particleX_i, cubicBSplineDiff_i, cubicBSplineDiff_j, e_x, mass_faceX_i);

//            A_X(cellIndex_i, cellIndex_j)+=Acomponent;
//          }
//        }
//      }

//      //Face Y
//      //Loop over particles in cell face Y for cellIndex
//      for (int particleIterator_i=0; particleIterator_i<noParticles_FaceY_i; particleIterator_i++)
//      {
//        for (int particleIterator_j=0; particleIterator_j<noParticles_FaceY_j; particleIterator_j++)
//        {
//          Particle* particleY_i=m_cellFacesY[cellIndex_i]->m_interpolationData[particleIterator_i]->m_particle;
//          Particle* particleY_j=m_cellFacesY[cellIndex_j]->m_interpolationData[particleIterator_j]->m_particle;

//          if (particleY_i==particleY_j)
//          {
//            //Calculate Aij component from p
//            //Get interpolation weights
//            Eigen::Vector3f cubicBSplineDiff_i=m_cellFacesY[cellIndex_i]->m_interpolationData[particleIterator_i]->m_cubicBSpline_Diff;
//            Eigen::Vector3f cubicBSplineDiff_j=m_cellFacesY[cellIndex_j]->m_interpolationData[particleIterator_j]->m_cubicBSpline_Diff;

//            float mass_faceY_i=m_cellFacesY[cellIndex_i]->m_mass;
//            float Acomponent=calcAComponent_DeviatoricVelocity(particleY_i, cubicBSplineDiff_i, cubicBSplineDiff_j, e_y, mass_faceY_i);

//            A_Y(cellIndex_i, cellIndex_j)+=Acomponent;
//          }
//        }
//      }

//      //Face Z
//      //Loop over particles in cell face Y for cellIndex
//      for (int particleIterator_i=0; particleIterator_i<noParticles_FaceZ_i; particleIterator_i++)
//      {
//        for (int particleIterator_j=0; particleIterator_j<noParticles_FaceZ_j; particleIterator_j++)
//        {
//          Particle* particleZ_i=m_cellFacesZ[cellIndex_i]->m_interpolationData[particleIterator_i]->m_particle;
//          Particle* particleZ_j=m_cellFacesZ[cellIndex_j]->m_interpolationData[particleIterator_j]->m_particle;

//          if (particleZ_i==particleZ_j)
//          {
//            //Calculate Aij component from p
//            //Get interpolation weights
//            Eigen::Vector3f cubicBSplineDiff_i=m_cellFacesZ[cellIndex_i]->m_interpolationData[particleIterator_i]->m_cubicBSpline_Diff;
//            Eigen::Vector3f cubicBSplineDiff_j=m_cellFacesZ[cellIndex_j]->m_interpolationData[particleIterator_j]->m_cubicBSpline_Diff;

//            float mass_faceZ_i=m_cellFacesZ[cellIndex_i]->m_mass;
//            float Acomponent=calcAComponent_DeviatoricVelocity(particleZ_i, cubicBSplineDiff_i, cubicBSplineDiff_j, e_z, mass_faceZ_i);

//            A_Z(cellIndex_i, cellIndex_j)+=Acomponent;
//          }
//        }
//      }
//    }
//  }

//  omp_set_nested(0);
  
  //Call MINRES to calculate new velocity values
  Eigen::VectorXf solution_X(m_totNoCells);
  solution_X.setZero();
  Eigen::VectorXf solution_Y(m_totNoCells);
  solution_Y.setZero();
  Eigen::VectorXf solution_Z(m_totNoCells);
  solution_Z.setZero();
  Eigen::MatrixXf emptyPreconditioner;
//  float shift=(-1.0);
  float shift=(0.0);
  float tolerance=0.0000001;
  int maxNoLoops=20;

  Eigen::MatrixXf A_X_trans=_A_X.transpose();
  Eigen::MatrixXf test=_A_X-A_X_trans;



  MathFunctions::MinRes(_A_X, _bVector_X, solution_X, emptyPreconditioner, shift, maxNoLoops, tolerance, false);
  MathFunctions::MinRes(_A_Y, _bVector_Y, solution_Y, emptyPreconditioner, shift, maxNoLoops, tolerance, false);
  MathFunctions::MinRes(_A_Z, _bVector_Z, solution_Z, emptyPreconditioner, shift, maxNoLoops, tolerance, false);


  //Read in solutions
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    m_cellFacesX[cellIndex]->m_velocity=solution_X[cellIndex];
    m_cellFacesY[cellIndex]->m_velocity=solution_Y[cellIndex];
    m_cellFacesZ[cellIndex]->m_velocity=solution_Z[cellIndex];
  }

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::searchCellsForCommonParticle(unsigned int _particleId, CellFace* _cellFace, unsigned int &o_particleIndexInFace, bool &o_isFound)
{
//  Particle* sameParticlePointer=nullptr;
  o_isFound=false;

  //Get number of particles in cell
  int noParticlesInList=_cellFace->m_interpolationData.size();

  //Check that id contained in list of j
  unsigned int particleId_Min=_cellFace->m_interpolationData[0]->m_particle->getId();

  if (_particleId>particleId_Min)
  {
    //Check if smaller than particle id max
    unsigned int particleId_Max=_cellFace->m_interpolationData[noParticlesInList-1]->m_particle->getId();

    if (_particleId<particleId_Max)
    {
      //Search through list
      int lowerBound=0;
      int upperBound=(noParticlesInList-1);
      int middleIndex=0;

      while (lowerBound<=upperBound)
      {
        middleIndex=lowerBound+((upperBound-lowerBound)/2);

        //Get particle id
        unsigned int midParticleId=_cellFace->m_interpolationData[middleIndex]->m_particle->getId();

        if (_particleId==midParticleId)
        {
//          o_particle=_cellFace->m_interpolationData[middleIndex]->m_particle;
          o_particleIndexInFace=middleIndex;
          o_isFound=true;
          break;
        }

        if (_particleId<midParticleId)
        {
          upperBound=middleIndex-1;
        }

        if (_particleId>midParticleId)
        {
          lowerBound=middleIndex+1;
        }
      }
    }

    if (_particleId==particleId_Max)
    {
//      o_particle=_cellFace->m_interpolationData[noParticlesInList-1]->m_particle;
      o_particleIndexInFace=noParticlesInList-1;
      o_isFound=true;
    }
  }

  if (_particleId==particleId_Min)
  {
//    o_particle=_cellFace->m_interpolationData[0]->m_particle;
    o_particleIndexInFace=0;
    o_isFound=true;
  }

//  *o_particle=*sameParticlePointer;
}

//----------------------------------------------------------------------------------------------------------------------
