#include "Grid.h"


//----------------------------------------------------------------------------------------------------------------------

void Grid::interpolateParticleToGrid(Emitter *_emitter, bool _isFirstStep)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Loop over particles in emitter

    Find position in grid

    Loop over all cells it could affect +-2?

      Calculate cubic B splines

      Add to mass for all cell centres and faces


  Loop over particles again

    Find position in grid

    Loop over all cells it could affect +-2?

      Calculate cubic B splines and their differentials

      Add to other variables to be updated

      If first step
        Add to particle density

//      Add to force

//      Add to B component

//      Loop over all possible cells for each cell

//        Add to A component



  ------------------------------------------------------------------------------------------------------
  */

  int noParticles=_emitter->m_noParticles;

  //Set e_{a(i)} vectors
  Eigen::Vector3f e_x(1.0, 0.0, 0.0);
  Eigen::Vector3f e_y(0.0, 1.0, 0.0);
  Eigen::Vector3f e_z(0.0, 0.0, 1.0);

  //To calc position of particle, need origin of grid corner, not centre of first grid cell.
  float halfCellSize=m_cellSize/2.0;
  Eigen::Vector3f gridEdgePosition=m_origin;
  gridEdgePosition(0)-=halfCellSize;
  gridEdgePosition(1)-=halfCellSize;
  gridEdgePosition(2)-=halfCellSize;

#pragma omp parallel for
  for (int particleItr=0; particleItr<noParticles; particleItr++)
  {
    Particle* particlePtr=_emitter->m_particles[particleItr];
    Eigen::Vector3f particlePosition=particlePtr->getPosition();
    Eigen::Vector3i particleIndex=MathFunctions::getParticleGridCell(particlePosition, m_cellSize, gridEdgePosition);

    //Get particle mass
    float particleMass=particlePtr->getMass();

    //Loop over i+-2, j+-2, k+-2 to get cells that particle will contribute to
    int iParticle=particleIndex(0);
    int jParticle=particleIndex(1);
    int kParticle=particleIndex(2);

    for (int kIndex=(kParticle-2); kIndex<(kParticle+4); kIndex++)
    {
      for (int jIndex=(jParticle-2); jIndex<(jParticle+4); jIndex++)
      {
        for (int iIndex=(iParticle-2); iIndex<(iParticle+4); iIndex++)
        {
          //Check that not in outer cells or outside grid
          if (iIndex>=0 && iIndex<=(m_noCells-1) && jIndex>=0 && jIndex<=(m_noCells-1) && kIndex>=0 && kIndex<=(m_noCells-1))
          {
            //Get cell index
            int cellIndex=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex, m_noCells);

            //Get cubic B Spline weight
            float weightCentre=0.0;
            float weightX=0.0;
            float weightY=0.0;
            float weightZ=0.0;
            calcWeight_cubicBSpline(particlePosition, iIndex, jIndex, kIndex,
                                    weightCentre, weightX, weightY, weightZ);

            //Add to number of particles contributing to faces and centres
            if (weightX!=0)
            {
              m_cellFacesX[cellIndex]->m_noParticlesContributing+=1;
              m_cellFacesX[cellIndex]->m_mass+=(weightX*particleMass);
            }
            if (weightY!=0)
            {
              m_cellFacesY[cellIndex]->m_noParticlesContributing+=1;
              m_cellFacesY[cellIndex]->m_mass+=(weightY*particleMass);
            }
            if (weightZ!=0)
            {
              m_cellFacesZ[cellIndex]->m_noParticlesContributing+=1;
              m_cellFacesZ[cellIndex]->m_mass+=(weightZ*particleMass);
            }
            if (weightCentre!=0)
            {
              m_cellCentres[cellIndex]->m_noParticlesContributing+=1;
              m_cellCentres[cellIndex]->m_mass+=(weightCentre*particleMass);
            }

          }
        }
      }
    }
  }

  //Calculate particle density if first step
  if (_isFirstStep==true)
  {
    #pragma omp parallel for
    for (int particleItr=0; particleItr<noParticles; particleItr++)
    {
      Particle* particlePtr=_emitter->m_particles[particleItr];
      Eigen::Vector3f particlePosition=particlePtr->getPosition();
      Eigen::Vector3i particleIndex=MathFunctions::getParticleGridCell(particlePosition, m_cellSize, gridEdgePosition);

      //Loop over i+-2, j+-2, k+-2 to get cells that particle will contribute to
      int iParticle=particleIndex(0);
      int jParticle=particleIndex(1);
      int kParticle=particleIndex(2);

      //Calc cell volume
      float cellVolume=pow(m_cellSize,3);

      for (int kIndex=(kParticle-2); kIndex<(kParticle+4); kIndex++)
      {
        for (int jIndex=(jParticle-2); jIndex<(jParticle+4); jIndex++)
        {
          for (int iIndex=(iParticle-2); iIndex<(iParticle+4); iIndex++)
          {
            //Check that not in outer cells or outside grid
            if (iIndex>=0 && iIndex<=(m_noCells-1) && jIndex>=0 && jIndex<=(m_noCells-1) && kIndex>=0 && kIndex<=(m_noCells-1))
            {
              //Get cell index
              int cellIndex=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex, m_noCells);

              //Get cubic B Spline weight
              float weightCentre=0.0;
              float weightX=0.0;
              float weightY=0.0;
              float weightZ=0.0;
              calcWeight_cubicBSpline(particlePosition, iIndex, jIndex, kIndex,
                                      weightCentre, weightX, weightY, weightZ);

              //Get mass of cell centre
              float mass=m_cellCentres[cellIndex]->m_mass;

              //Add to density of particle
              float densityContrib=(mass*weightCentre)/cellVolume;

              particlePtr->addParticleDensity(densityContrib);

            }
          }
        }
      }

      //When all density contributions have been added, calc initial volume
      particlePtr->calcInitialVolume();
    }
  }


  //Loop over particles to rasterise particle data to grid
//  #pragma omp parallel for
  for (int particleItr=0; particleItr<noParticles; particleItr++)
  {
    //Get particle pointer
    Particle* particlePtr=_emitter->m_particles[particleItr];

    //Get particle variables
    float particleMass=0.0;
    Eigen::Vector3f particleVelocity;
    float particleHeatConductivity;
    float particleDetDeformGrad;
    float particleDetDeformGradElastic;
    float particleHeatCapacity;
    Phase particlePhase;
    float particleTemperature;
    float particleLameLambdaInv;

    particlePtr->getParticleData_CellFace(particleMass, particleVelocity, particlePhase);
    particlePtr->getParticleData_CellCentre(particleMass, particleDetDeformGrad, particleDetDeformGradElastic, particlePhase, particleTemperature, particleLameLambdaInv);

    if (particlePhase==Phase::Liquid)
    {
      particleHeatCapacity=_emitter->m_heatCapacityFluid;
      particleHeatConductivity=_emitter->m_heatConductivityFluid;
    }
    else
    {
      particleHeatCapacity=_emitter->m_heatCapacitySolid;
      particleHeatConductivity=_emitter->m_heatConductivitySolid;
    }

    //Get particle position in grid
    Eigen::Vector3f particlePosition=particlePtr->getPosition();
    Eigen::Vector3i particleIndex=MathFunctions::getParticleGridCell(particlePosition, m_cellSize, gridEdgePosition);

    //Loop over i+-2, j+-2, k+-2 to get cells that particle will contribute to
    int iParticle=particleIndex(0);
    int jParticle=particleIndex(1);
    int kParticle=particleIndex(2);

//    omp_set_nested(1);
#pragma omp parallel for
    for (int kIndex=(kParticle-2); kIndex<(kParticle+4); kIndex++)
    {
      for (int jIndex=(jParticle-2); jIndex<(jParticle+4); jIndex++)
      {
        for (int iIndex=(iParticle-2); iIndex<(iParticle+4); iIndex++)
        {
          //Check that not in outer cells or outside grid
          if (iIndex>=0 && iIndex<=(m_noCells-1) && jIndex>=0 && jIndex<=(m_noCells-1) && kIndex>=0 && kIndex<=(m_noCells-1))
          {
            //Get cell index
            int cellIndex=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex, m_noCells);

            //Get cubic B Spline weight
            float weightCentre=0.0;
            float weightX=0.0;
            float weightY=0.0;
            float weightZ=0.0;
            calcWeight_cubicBSpline(particlePosition, iIndex, jIndex, kIndex,
                                    weightCentre, weightX, weightY, weightZ);

            //Mass of cell faces and centre
            float mass_Centre=m_cellCentres[cellIndex]->m_mass;
            float mass_FaceX=m_cellFacesX[cellIndex]->m_mass;
            float mass_FaceY=m_cellFacesY[cellIndex]->m_mass;
            float mass_FaceZ=m_cellFacesZ[cellIndex]->m_mass;

            //Get number of particles contributing to cell centre and faces
            int noParticles_Centre=m_cellCentres[cellIndex]->m_noParticlesContributing;
            int noParticles_FaceX=m_cellFacesX[cellIndex]->m_noParticlesContributing;
            int noParticles_FaceY=m_cellFacesY[cellIndex]->m_noParticlesContributing;
            int noParticles_FaceZ=m_cellFacesZ[cellIndex]->m_noParticlesContributing;

            if (noParticles_FaceX>0)
            {
              m_cellFacesX[cellIndex]->m_velocity+=(particleMass*weightX*particleVelocity(0))/mass_FaceX;
              m_cellFacesX[cellIndex]->m_heatConductivity+=(particleMass*weightX*particleHeatConductivity)/mass_FaceX;
            }

            if (noParticles_FaceY>0)
            {
              m_cellFacesY[cellIndex]->m_velocity+=(particleMass*weightY*particleVelocity(1))/mass_FaceY;
              m_cellFacesY[cellIndex]->m_heatConductivity+=(particleMass*weightY*particleHeatConductivity)/mass_FaceY;
            }

            if (noParticles_FaceZ>0)
            {
              m_cellFacesZ[cellIndex]->m_velocity+=(particleMass*weightZ*particleVelocity(2))/mass_FaceZ;
              m_cellFacesZ[cellIndex]->m_heatConductivity+=(particleMass*weightZ*particleHeatConductivity)/mass_FaceZ;
            }

            if (noParticles_Centre>0)
            {
              //Centre variables
              m_cellCentres[cellIndex]->m_detDeformationGrad+=(particleMass*weightCentre*particleDetDeformGrad)/mass_Centre;
              m_cellCentres[cellIndex]->m_detDeformationGradElastic+=(particleMass*weightCentre*particleDetDeformGradElastic)/mass_Centre;
              m_cellCentres[cellIndex]->m_heatCapacity+=(particleMass*weightCentre*particleHeatCapacity)/mass_Centre;
              m_cellCentres[cellIndex]->m_temperature+=(particleMass*weightCentre*particleTemperature)/mass_Centre;
              m_cellCentres[cellIndex]->m_lameLambdaInverse+=(particleMass*weightCentre*particleLameLambdaInv)/mass_Centre;
            }


            //Calculate force, B component and A row components for deviatoric velocity
            if (noParticles_FaceX>0 || noParticles_FaceY>0 || noParticles_FaceZ>0)
            {
              calcDeviatoricContributions(particlePtr, cellIndex, iIndex, jIndex, kIndex, weightX, weightY, weightZ);
            }

          }
        }
      }
    }

  }

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcWeight_cubicBSpline(Eigen::Vector3f _particlePosition, int _iIndex, int _jIndex, int _kIndex,
                                   float &o_weightCentre, float &o_weightFaceX, float &o_weightFaceY, float &o_weightFaceZ)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Calculate interpolation weights for the particle and the cell with given _i, _j,_k

  Calculate i*, ix, iy, iz position vectors for the cell in question

  Calculate x=xp-xia values for these

  Pass in x and get interpolation weights:
       cubicBspline

  Return results
  ------------------------------------------------------------------------------------------------------
  */

  //Cell position
  float xPos=(_iIndex*m_cellSize)+m_origin(0);
  float yPos=(_jIndex*m_cellSize)+m_origin(1);
  float zPos=(_kIndex*m_cellSize)+m_origin(2);

  //Position vectors for centre and faces
  Eigen::Vector3f centreVector(xPos, yPos, zPos);

  float halfCellSize=m_cellSize/2.0;
  Eigen::Vector3f faceXVector(xPos-halfCellSize, yPos, zPos);
  Eigen::Vector3f faceYVector(xPos, yPos-halfCellSize, zPos);
  Eigen::Vector3f faceZVector(xPos, yPos, zPos-halfCellSize);


  //Calculate posDifference for each face and cell centre
  Eigen::Vector3f centrePosDiff=_particlePosition-centreVector;
  Eigen::Vector3f faceXPosDiff=_particlePosition-faceXVector;
  Eigen::Vector3f faceYPosDiff=_particlePosition-faceYVector;
  Eigen::Vector3f faceZPosDiff=_particlePosition-faceZVector;

  //Centre
  float NxCentre_cubicBS=MathFunctions::calcCubicBSpline(centrePosDiff(0)/m_cellSize);
  float NyCentre_cubicBS=MathFunctions::calcCubicBSpline(centrePosDiff(1)/m_cellSize);
  float NzCentre_cubicBS=MathFunctions::calcCubicBSpline(centrePosDiff(2)/m_cellSize);

  //FaceX
  float NxFaceX_cubicBS=MathFunctions::calcCubicBSpline(faceXPosDiff(0)/m_cellSize);
  float NyFaceX_cubicBS=MathFunctions::calcCubicBSpline(faceXPosDiff(1)/m_cellSize);
  float NzFaceX_cubicBS=MathFunctions::calcCubicBSpline(faceXPosDiff(2)/m_cellSize);

  //FaceY
  float NxFaceY_cubicBS=MathFunctions::calcCubicBSpline(faceYPosDiff(0)/m_cellSize);
  float NyFaceY_cubicBS=MathFunctions::calcCubicBSpline(faceYPosDiff(1)/m_cellSize);
  float NzFaceY_cubicBS=MathFunctions::calcCubicBSpline(faceYPosDiff(2)/m_cellSize);

  //FaceZ
  float NxFaceZ_cubicBS=MathFunctions::calcCubicBSpline(faceZPosDiff(0)/m_cellSize);
  float NyFaceZ_cubicBS=MathFunctions::calcCubicBSpline(faceZPosDiff(1)/m_cellSize);
  float NzFaceZ_cubicBS=MathFunctions::calcCubicBSpline(faceZPosDiff(2)/m_cellSize);

  //Return interpolation weights
  o_weightCentre=NxCentre_cubicBS*NyCentre_cubicBS*NzCentre_cubicBS;
  o_weightFaceX=NxFaceX_cubicBS*NyFaceX_cubicBS*NzFaceX_cubicBS;
  o_weightFaceY=NxFaceY_cubicBS*NyFaceY_cubicBS*NzFaceY_cubicBS;
  o_weightFaceZ=NxFaceZ_cubicBS*NyFaceZ_cubicBS*NzFaceZ_cubicBS;



}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcWeight_cubicBSpline_Diff(Eigen::Vector3f _particlePosition, int _iIndex, int _jIndex, int _kIndex,
                                        Eigen::Vector3f &o_weightFaceX, Eigen::Vector3f &o_weightFaceY, Eigen::Vector3f &o_weightFaceZ)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Calculate interpolation weights for the particle and the cell with given _i, _j,_k

  Calculate i*, ix, iy, iz position vectors for the cell in question

  Calculate x=xp-xia values for these

  Pass in x and get interpolation weights:
       cubicBspline

  Return results
  ------------------------------------------------------------------------------------------------------
  */

  //Cell position
  float xPos=(_iIndex*m_cellSize)+m_origin(0);
  float yPos=(_jIndex*m_cellSize)+m_origin(1);
  float zPos=(_kIndex*m_cellSize)+m_origin(2);

  //Position vectors for centre and faces
  Eigen::Vector3f centreVector(xPos, yPos, zPos);

  float halfCellSize=m_cellSize/2.0;
  Eigen::Vector3f faceXVector(xPos-halfCellSize, yPos, zPos);
  Eigen::Vector3f faceYVector(xPos, yPos-halfCellSize, zPos);
  Eigen::Vector3f faceZVector(xPos, yPos, zPos-halfCellSize);

  //Calculate posDifference for each face and cell centre
  Eigen::Vector3f centrePosDiff=_particlePosition-centreVector;
  Eigen::Vector3f faceXPosDiff=_particlePosition-faceXVector;
  Eigen::Vector3f faceYPosDiff=_particlePosition-faceYVector;
  Eigen::Vector3f faceZPosDiff=_particlePosition-faceZVector;

  //Need to get the non-differentiated values as well
  //Centre
  float NxCentre_cubicBS=MathFunctions::calcCubicBSpline(centrePosDiff(0)/m_cellSize);
  float NyCentre_cubicBS=MathFunctions::calcCubicBSpline(centrePosDiff(1)/m_cellSize);
  float NzCentre_cubicBS=MathFunctions::calcCubicBSpline(centrePosDiff(2)/m_cellSize);

  //Face X
  float dNx_cubicBS_FaceX=MathFunctions::calcCubicBSpline_Diff(faceXPosDiff(0)/m_cellSize);
  float dNy_cubicBS_FaceX=MathFunctions::calcCubicBSpline_Diff(faceXPosDiff(1)/m_cellSize);
  float dNz_cubicBS_FaceX=MathFunctions::calcCubicBSpline_Diff(faceXPosDiff(2)/m_cellSize);

  o_weightFaceX(0)=dNx_cubicBS_FaceX*NyCentre_cubicBS*NzCentre_cubicBS;
  o_weightFaceX(1)=dNy_cubicBS_FaceX*NxCentre_cubicBS*NzCentre_cubicBS;
  o_weightFaceX(2)=dNz_cubicBS_FaceX*NxCentre_cubicBS*NyCentre_cubicBS;

  o_weightFaceX*=(1.0/m_cellSize); ///Not sure about this part?
//  o_weightFaceX*=(-1.0/m_cellSize); ///Not sure about this part?

  //Face Y
  float dNx_cubicBS_FaceY=MathFunctions::calcCubicBSpline_Diff(faceYPosDiff(0)/m_cellSize);
  float dNy_cubicBS_FaceY=MathFunctions::calcCubicBSpline_Diff(faceYPosDiff(1)/m_cellSize);
  float dNz_cubicBS_FaceY=MathFunctions::calcCubicBSpline_Diff(faceYPosDiff(2)/m_cellSize);

  o_weightFaceY(0)=dNx_cubicBS_FaceY*NyCentre_cubicBS*NzCentre_cubicBS;
  o_weightFaceY(1)=dNy_cubicBS_FaceY*NxCentre_cubicBS*NzCentre_cubicBS;
  o_weightFaceY(2)=dNz_cubicBS_FaceY*NxCentre_cubicBS*NyCentre_cubicBS;

  o_weightFaceY*=(1.0/m_cellSize); ///Not sure about this part?
//  o_weightFaceY*=(-1.0/m_cellSize); ///Not sure about this part?

  //Face Z
  float dNx_cubicBS_FaceZ=MathFunctions::calcCubicBSpline_Diff(faceZPosDiff(0)/m_cellSize);
  float dNy_cubicBS_FaceZ=MathFunctions::calcCubicBSpline_Diff(faceZPosDiff(1)/m_cellSize);
  float dNz_cubicBS_FaceZ=MathFunctions::calcCubicBSpline_Diff(faceZPosDiff(2)/m_cellSize);

  o_weightFaceZ(0)=dNx_cubicBS_FaceZ*NyCentre_cubicBS*NzCentre_cubicBS;
  o_weightFaceZ(1)=dNy_cubicBS_FaceZ*NxCentre_cubicBS*NzCentre_cubicBS;
  o_weightFaceZ(2)=dNz_cubicBS_FaceZ*NxCentre_cubicBS*NyCentre_cubicBS;

  o_weightFaceZ*=(1.0/m_cellSize); ///Not sure about this part?
//  o_weightFaceZ*=(-1.0/m_cellSize); ///Not sure about this part?

}

void Grid::classifyCells_New()
{
  /* Outline - Current setup might be time consuming since two loops. Could rectify this for bounding box, but not for level sets I think
  ----------------------------------------------------------------------------------------------------------------------
  Loop over all cell faces
    Check faces against collision
      For now use bounding box, so colliding if i<2||>n-2, j<2||>n-2, k<2||>n-2

  Loop over all cell centres
    Check if 3 faces are colliding - If cell centre is i<n-1, j<n-1 or k<n-1 then check the faces of the nearest neighbour
    cells in each of the directions as well.
      If all colliding - colliding
      If not all and no particles - empty
      Otherwise - interior
    Set heat source temperature for colliding cells with kIndex==0. Ie. heat source element is the k=0 plane.
    Set ambient temperature to empty cells

  ----------------------------------------------------------------------------------------------------------------------
  */

  //Loop over cell faces - This loop could be made smaller when just checking the outer cells.
  //But this is possibly easier to thread
#pragma omp parallel for
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    //Test parallel
//    printf("The parallel region is executed by thread %d\n", omp_get_thread_num());

    //Find cell index. Will be same for the other faces
    int iIndex=m_cellFacesX[cellIndex]->m_iIndex;
    int jIndex=m_cellFacesX[cellIndex]->m_jIndex;
    int kIndex=m_cellFacesX[cellIndex]->m_kIndex;

    //This checks whether cell faces belong to outer cells. To set collision cells to be the outer rim of cells
    if (iIndex==0 || iIndex==(m_noCells-1) || jIndex==0 || jIndex==(m_noCells-1) || kIndex==0 || kIndex==(m_noCells-1) )
    {
      //Set faces to colliding
      m_cellFacesX[cellIndex]->m_state=State::Colliding;
      m_cellFacesY[cellIndex]->m_state=State::Colliding;
      m_cellFacesZ[cellIndex]->m_state=State::Colliding;
    }

    //Also need to set cell faces adjacent to the outer cells to colliding. This must be done separately for each cell
    if (iIndex==1)
    {
      m_cellFacesX[cellIndex]->m_state=State::Colliding;
    }
    if (jIndex==1)
    {
      m_cellFacesY[cellIndex]->m_state=State::Colliding;
    }
    if (kIndex==1)
    {
      m_cellFacesZ[cellIndex]->m_state=State::Colliding;
    }
  }

  //This step will work for level set collisions as well.
  //Loop over all cells again to check which cell centres are collding
  //Seems inefficient.
#pragma omp parallel for
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    //Test parallel
//    printf("The parallel region is executed by thread %d\n", omp_get_thread_num());

    //Find cell index.
    int iIndex=m_cellCentres[cellIndex]->m_iIndex;
    int jIndex=m_cellCentres[cellIndex]->m_jIndex;
    int kIndex=m_cellCentres[cellIndex]->m_kIndex;

    //Get indices of faces in the positive ijk directions
    int cellIndex_i1jk=MathFunctions::getVectorIndex(iIndex+1, jIndex, kIndex, m_noCells);
    int cellIndex_ij1k=MathFunctions::getVectorIndex(iIndex, jIndex+1, kIndex, m_noCells);
    int cellIndex_ijk1=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex+1, m_noCells);

    //Face X
    //Check if lower x face colliding
    if (m_cellFacesX[cellIndex]->m_state!=State::Colliding)
    {
      //Check whether empty or not
      int noParticlesInCellCentre=m_cellCentres[cellIndex]->m_noParticlesContributing;
      int noParticlesInCellFaceX_1=m_cellFacesX[cellIndex]->m_noParticlesContributing;
      int noParticlesInCellFaceY_1=m_cellFacesY[cellIndex]->m_noParticlesContributing;
      int noParticlesInCellFaceZ_1=m_cellFacesZ[cellIndex]->m_noParticlesContributing;

      //Get particle number for upper faces of cell, unless outermost cells in grid
      int noParticlesInCellFaceX1=0;
      int noParticlesInCellFaceY1=0;
      int noParticlesInCellFaceZ1=0;

      if (iIndex!=(m_noCells-1))
      {
        noParticlesInCellFaceX1=m_cellFacesX[cellIndex_i1jk]->m_noParticlesContributing;
      }
      if (jIndex!=(m_noCells-1))
      {
        noParticlesInCellFaceY1=m_cellFacesY[cellIndex_ij1k]->m_noParticlesContributing;
      }
      if (kIndex!=(m_noCells-1))
      {
        noParticlesInCellFaceZ1=m_cellFacesZ[cellIndex_ijk1]->m_noParticlesContributing;
      }

      //If cell centre and all faces belonging to cells have particles affecting it, then cell is interior
      if (noParticlesInCellCentre>m_noParticlesThreshold
          && noParticlesInCellFaceX1>m_noParticlesThreshold
          && noParticlesInCellFaceX_1>m_noParticlesThreshold
          && noParticlesInCellFaceY1>m_noParticlesThreshold
          && noParticlesInCellFaceY_1>m_noParticlesThreshold
          && noParticlesInCellFaceZ1>m_noParticlesThreshold
          && noParticlesInCellFaceZ_1>m_noParticlesThreshold)
      {
        m_cellCentres[cellIndex]->m_state=State::Interior;
      }
      //Otherwise the cell is empty
      else
      {
        m_cellCentres[cellIndex]->m_state=State::Empty;
        m_cellCentres[cellIndex]->m_temperature=m_ambientTemperature;
      }

      //Set face to empty as well
      if (noParticlesInCellFaceX_1<=m_noParticlesThreshold)
      {
        m_cellFacesX[cellIndex]->m_state=State::Empty;
      }

      //Go to next cellIndex
      continue;
    }

    //Check upper x face as well
    if (iIndex<(m_noCells-1))
    {
      if (m_cellFacesX[cellIndex_i1jk]->m_state!=State::Colliding)
      {
        //Check whether empty or not
        int noParticlesInCellCentre=m_cellCentres[cellIndex]->m_noParticlesContributing;
        int noParticlesInCellFaceX_1=m_cellFacesX[cellIndex]->m_noParticlesContributing;
        int noParticlesInCellFaceY_1=m_cellFacesY[cellIndex]->m_noParticlesContributing;
        int noParticlesInCellFaceZ_1=m_cellFacesZ[cellIndex]->m_noParticlesContributing;

        //Get particle number for upper faces of cell, unless outermost cells in grid
        int noParticlesInCellFaceX1=0;
        int noParticlesInCellFaceY1=0;
        int noParticlesInCellFaceZ1=0;

        if (iIndex!=(m_noCells-1))
        {
          noParticlesInCellFaceX1=m_cellFacesX[cellIndex_i1jk]->m_noParticlesContributing;
        }
        if (jIndex!=(m_noCells-1))
        {
          noParticlesInCellFaceY1=m_cellFacesY[cellIndex_ij1k]->m_noParticlesContributing;
        }
        if (kIndex!=(m_noCells-1))
        {
          noParticlesInCellFaceZ1=m_cellFacesZ[cellIndex_ijk1]->m_noParticlesContributing;
        }

        //If cell centre and all faces belonging to cells have particles affecting it, then cell is interior

        if (noParticlesInCellCentre>m_noParticlesThreshold
            && noParticlesInCellFaceX1>m_noParticlesThreshold
            && noParticlesInCellFaceX_1>m_noParticlesThreshold
            && noParticlesInCellFaceY1>m_noParticlesThreshold
            && noParticlesInCellFaceY_1>m_noParticlesThreshold
            && noParticlesInCellFaceZ1>m_noParticlesThreshold
            && noParticlesInCellFaceZ_1>m_noParticlesThreshold)
        {
          m_cellCentres[cellIndex]->m_state=State::Interior;
        }
        //Otherwise the cell is empty
        else
        {
          m_cellCentres[cellIndex]->m_state=State::Empty;
          m_cellCentres[cellIndex]->m_temperature=m_ambientTemperature;
        }

        //Go to next cellIndex
        continue;
      }
    }

    //Face Y
    //Check if lower x face colliding
    if (m_cellFacesY[cellIndex]->m_state!=State::Colliding)
    {
      //Check whether empty or not
      int noParticlesInCellCentre=m_cellCentres[cellIndex]->m_noParticlesContributing;
      int noParticlesInCellFaceX_1=m_cellFacesX[cellIndex]->m_noParticlesContributing;
      int noParticlesInCellFaceY_1=m_cellFacesY[cellIndex]->m_noParticlesContributing;
      int noParticlesInCellFaceZ_1=m_cellFacesZ[cellIndex]->m_noParticlesContributing;

      //Get particle number for upper faces of cell, unless outermost cells in grid
      int noParticlesInCellFaceX1=0;
      int noParticlesInCellFaceY1=0;
      int noParticlesInCellFaceZ1=0;

      if (iIndex!=(m_noCells-1))
      {
        noParticlesInCellFaceX1=m_cellFacesX[cellIndex_i1jk]->m_noParticlesContributing;
      }
      if (jIndex!=(m_noCells-1))
      {
        noParticlesInCellFaceY1=m_cellFacesY[cellIndex_ij1k]->m_noParticlesContributing;
      }
      if (kIndex!=(m_noCells-1))
      {
        noParticlesInCellFaceZ1=m_cellFacesZ[cellIndex_ijk1]->m_noParticlesContributing;
      }

      //If cell centre and all faces belonging to cells have particles affecting it, then cell is interior
      if (noParticlesInCellCentre>m_noParticlesThreshold
          && noParticlesInCellFaceX1>m_noParticlesThreshold
          && noParticlesInCellFaceX_1>m_noParticlesThreshold
          && noParticlesInCellFaceY1>m_noParticlesThreshold
          && noParticlesInCellFaceY_1>m_noParticlesThreshold
          && noParticlesInCellFaceZ1>m_noParticlesThreshold
          && noParticlesInCellFaceZ_1>m_noParticlesThreshold)
      {
        m_cellCentres[cellIndex]->m_state=State::Interior;
      }
      //Otherwise the cell is empty
      else
      {
        m_cellCentres[cellIndex]->m_state=State::Empty;
        m_cellCentres[cellIndex]->m_temperature=m_ambientTemperature;
      }

      //Set face to empty as well
      if (noParticlesInCellFaceY_1<=m_noParticlesThreshold)
      {
        m_cellFacesY[cellIndex]->m_state=State::Empty;
      }

      //Go to next cellIndex
      continue;
    }

    //Check upper y face as well
    if (jIndex<(m_noCells-1))
    {
      if (m_cellFacesY[cellIndex_ij1k]->m_state!=State::Colliding)
      {
        //Check whether empty or not
        int noParticlesInCellCentre=m_cellCentres[cellIndex]->m_noParticlesContributing;
        int noParticlesInCellFaceX_1=m_cellFacesX[cellIndex]->m_noParticlesContributing;
        int noParticlesInCellFaceY_1=m_cellFacesY[cellIndex]->m_noParticlesContributing;
        int noParticlesInCellFaceZ_1=m_cellFacesZ[cellIndex]->m_noParticlesContributing;

        //Get particle number for upper faces of cell, unless outermost cells in grid
        int noParticlesInCellFaceX1=0;
        int noParticlesInCellFaceY1=0;
        int noParticlesInCellFaceZ1=0;

        if (iIndex!=(m_noCells-1))
        {
          noParticlesInCellFaceX1=m_cellFacesX[cellIndex_i1jk]->m_noParticlesContributing;
        }
        if (jIndex!=(m_noCells-1))
        {
          noParticlesInCellFaceY1=m_cellFacesY[cellIndex_ij1k]->m_noParticlesContributing;
        }
        if (kIndex!=(m_noCells-1))
        {
          noParticlesInCellFaceZ1=m_cellFacesZ[cellIndex_ijk1]->m_noParticlesContributing;
        }

        //If cell centre and all faces belonging to cells have particles affecting it, then cell is interior
        if (noParticlesInCellCentre>m_noParticlesThreshold
            && noParticlesInCellFaceX1>m_noParticlesThreshold
            && noParticlesInCellFaceX_1>m_noParticlesThreshold
            && noParticlesInCellFaceY1>m_noParticlesThreshold
            && noParticlesInCellFaceY_1>m_noParticlesThreshold
            && noParticlesInCellFaceZ1>m_noParticlesThreshold
            && noParticlesInCellFaceZ_1>m_noParticlesThreshold)
        {
          m_cellCentres[cellIndex]->m_state=State::Interior;
        }
        //Otherwise the cell is empty
        else
        {
          m_cellCentres[cellIndex]->m_state=State::Empty;
          m_cellCentres[cellIndex]->m_temperature=m_ambientTemperature;
        }

        //Go to next cellIndex
        continue;
      }
    }

    //Face Z
    //Check if lower x face colliding
    if (m_cellFacesZ[cellIndex]->m_state!=State::Colliding)
    {
      //Check whether empty or not
      int noParticlesInCellCentre=m_cellCentres[cellIndex]->m_noParticlesContributing;
      int noParticlesInCellFaceX_1=m_cellFacesX[cellIndex]->m_noParticlesContributing;
      int noParticlesInCellFaceY_1=m_cellFacesY[cellIndex]->m_noParticlesContributing;
      int noParticlesInCellFaceZ_1=m_cellFacesZ[cellIndex]->m_noParticlesContributing;

      //Get particle number for upper faces of cell, unless outermost cells in grid
      int noParticlesInCellFaceX1=0;
      int noParticlesInCellFaceY1=0;
      int noParticlesInCellFaceZ1=0;

      if (iIndex!=(m_noCells-1))
      {
        noParticlesInCellFaceX1=m_cellFacesX[cellIndex_i1jk]->m_noParticlesContributing;
      }
      if (jIndex!=(m_noCells-1))
      {
        noParticlesInCellFaceY1=m_cellFacesY[cellIndex_ij1k]->m_noParticlesContributing;
      }
      if (kIndex!=(m_noCells-1))
      {
        noParticlesInCellFaceZ1=m_cellFacesZ[cellIndex_ijk1]->m_noParticlesContributing;
      }

      //If cell centre and all faces belonging to cells have particles affecting it, then cell is interior
      if (noParticlesInCellCentre>m_noParticlesThreshold
          && noParticlesInCellFaceX1>m_noParticlesThreshold
          && noParticlesInCellFaceX_1>m_noParticlesThreshold
          && noParticlesInCellFaceY1>m_noParticlesThreshold
          && noParticlesInCellFaceY_1>m_noParticlesThreshold
          && noParticlesInCellFaceZ1>m_noParticlesThreshold
          && noParticlesInCellFaceZ_1>m_noParticlesThreshold)
      {
        m_cellCentres[cellIndex]->m_state=State::Interior;
      }
      //Otherwise the cell is empty
      else
      {
        m_cellCentres[cellIndex]->m_state=State::Empty;
        m_cellCentres[cellIndex]->m_temperature=m_ambientTemperature;
      }

      //Set face to empty as well
      if (noParticlesInCellFaceZ_1<=m_noParticlesThreshold)
      {
        m_cellFacesZ[cellIndex]->m_state=State::Empty;
      }

      //Go to next cellIndex
      continue;
    }

    //Check upper x face as well
    if (kIndex<(m_noCells-1))
    {
      if (m_cellFacesZ[cellIndex_ijk1]->m_state!=State::Colliding)
      {
        //Check whether empty or not
        int noParticlesInCellCentre=m_cellCentres[cellIndex]->m_noParticlesContributing;
        int noParticlesInCellFaceX_1=m_cellFacesX[cellIndex]->m_noParticlesContributing;
        int noParticlesInCellFaceY_1=m_cellFacesY[cellIndex]->m_noParticlesContributing;
        int noParticlesInCellFaceZ_1=m_cellFacesZ[cellIndex]->m_noParticlesContributing;

        //Get particle number for upper faces of cell, unless outermost cells in grid
        int noParticlesInCellFaceX1=0;
        int noParticlesInCellFaceY1=0;
        int noParticlesInCellFaceZ1=0;

        if (iIndex!=(m_noCells-1))
        {
          noParticlesInCellFaceX1=m_cellFacesX[cellIndex_i1jk]->m_noParticlesContributing;
        }
        if (jIndex!=(m_noCells-1))
        {
          noParticlesInCellFaceY1=m_cellFacesY[cellIndex_ij1k]->m_noParticlesContributing;
        }
        if (kIndex!=(m_noCells-1))
        {
          noParticlesInCellFaceZ1=m_cellFacesZ[cellIndex_ijk1]->m_noParticlesContributing;
        }

        //If cell centre and all faces belonging to cells have particles affecting it, then cell is interior
        if (noParticlesInCellCentre>m_noParticlesThreshold
            && noParticlesInCellFaceX1>m_noParticlesThreshold
            && noParticlesInCellFaceX_1>m_noParticlesThreshold
            && noParticlesInCellFaceY1>m_noParticlesThreshold
            && noParticlesInCellFaceY_1>m_noParticlesThreshold
            && noParticlesInCellFaceZ1>m_noParticlesThreshold
            && noParticlesInCellFaceZ_1>m_noParticlesThreshold)
        {
          m_cellCentres[cellIndex]->m_state=State::Interior;
        }
        //Otherwise the cell is empty
        else
        {
          m_cellCentres[cellIndex]->m_state=State::Empty;
          m_cellCentres[cellIndex]->m_temperature=m_ambientTemperature;
        }

        //Go to next cellIndex
        continue;
      }
    }

    //This section will not be reached if faces that are non-colliding are found
    //Set temperatures for colliding cells that are colliding with a heat source object
    //Heat source object set to j=0 plane, ie. jIndex==0
    if (jIndex==0)
    {
      m_cellCentres[cellIndex]->m_temperature=m_heatSourceTemperature;
    }
    //Need to set empty collision cells to ambient temperature
    else if (m_cellCentres[cellIndex]->m_noParticlesContributing==0)
    {
      m_cellCentres[cellIndex]->m_temperature=m_ambientTemperature;
    }

  }
}
