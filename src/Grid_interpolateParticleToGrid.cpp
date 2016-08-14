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

  int noParticles=_emitter->getNoParticles();

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
            }
            if (weightY!=0)
            {
              m_cellFacesY[cellIndex]->m_noParticlesContributing+=1;
            }
            if (weightZ!=0)
            {
              m_cellFacesZ[cellIndex]->m_noParticlesContributing+=1;
            }
            if (weightCentre!=0)
            {
              m_cellCentres[cellIndex]->m_noParticlesContributing+=1;
            }

            //Add to mass
            m_cellCentres[cellIndex]->m_mass+=(weightCentre*particleMass);
            m_cellFacesX[cellIndex]->m_mass+=(weightX*particleMass);
            m_cellFacesY[cellIndex]->m_mass+=(weightY*particleMass);
            m_cellFacesZ[cellIndex]->m_mass+=(weightZ*particleMass);

          }
        }
      }
    }
  }

  //Calculate particle density if first step
  if (_isFirstStep==true)
  {
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

              //Get mass of faces
              float massX=m_cellFacesX[cellIndex]->m_mass;
              float massY=m_cellFacesY[cellIndex]->m_mass;
              float massZ=m_cellFacesZ[cellIndex]->m_mass;

              //Add to density of particle
              float densityContrib=0.0;
              densityContrib+=((massX*weightX)/cellVolume);
              densityContrib+=((massY*weightY)/cellVolume);
              densityContrib+=((massZ*weightZ)/cellVolume);

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

            //Face variables
            m_cellFacesX[cellIndex]->m_velocity+=(particleMass*weightX*particleVelocity(0))/mass_FaceX;
            m_cellFacesX[cellIndex]->m_heatConductivity+=(particleMass*weightX*particleHeatConductivity)/mass_FaceX;
            m_cellFacesY[cellIndex]->m_velocity+=(particleMass*weightY*particleVelocity(1))/mass_FaceY;
            m_cellFacesY[cellIndex]->m_heatConductivity+=(particleMass*weightY*particleHeatConductivity)/mass_FaceY;
            m_cellFacesZ[cellIndex]->m_velocity+=(particleMass*weightZ*particleVelocity(2))/mass_FaceZ;
            m_cellFacesZ[cellIndex]->m_heatConductivity+=(particleMass*weightZ*particleHeatConductivity)/mass_FaceZ;

            //Centre variables
            m_cellCentres[cellIndex]->m_detDeformationGrad+=(particleMass*weightCentre*particleDetDeformGrad)/mass_Centre;
            m_cellCentres[cellIndex]->m_detDeformationGradElastic+=(particleMass*weightCentre*particleDetDeformGradElastic)/mass_Centre;
            m_cellCentres[cellIndex]->m_heatCapacity+=(particleMass*weightCentre*particleHeatCapacity)/mass_Centre;
            m_cellCentres[cellIndex]->m_temperature+=(particleMass*weightCentre*particleTemperature)/mass_Centre;
            m_cellCentres[cellIndex]->m_lameLambdaInverse+=(particleMass*weightCentre*particleLameLambdaInv)/mass_Centre;


            //Calculate force, B component and A row components for deviatoric velocity
            calcDeviatoricContributions(particlePtr, cellIndex, iIndex, jIndex, kIndex, weightX, weightY, weightZ);

          }
        }
      }
    }

  }

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcDeviatoricContributions(Particle *_particle, int _cellIndex, int _iIndex, int _jIndex, int _kIndex,
                                       float _weightX, float _weightY, float _weightZ)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Calc differentiated weight

  Add to force

  Add to B component

  Loop over all possible cells for each cell

    Add to A component
  ------------------------------------------------------------------------------------------------------
  */

  //Set e_{a(i)} vectors
  Eigen::Vector3f e_x(1.0, 0.0, 0.0);
  Eigen::Vector3f e_y(0.0, 1.0, 0.0);
  Eigen::Vector3f e_z(0.0, 0.0, 1.0);

  Eigen::Vector3f weightDiff_FaceX;
  Eigen::Vector3f weightDiff_FaceY;
  Eigen::Vector3f weightDiff_FaceZ;

  //Get particle position
  Eigen::Vector3f particlePosition=_particle->getPosition();

  //Get differentiated weights
  calcWeight_cubicBSpline_Diff(particlePosition, _iIndex, _jIndex, _kIndex, weightDiff_FaceX, weightDiff_FaceY, weightDiff_FaceZ);

  //Get deviatoric forces
  float force_FaceX=calcDeviatoricForce_New(_particle, e_x, weightDiff_FaceX);
  float force_FaceY=calcDeviatoricForce_New(_particle, e_y, weightDiff_FaceY);
  float force_FaceZ=calcDeviatoricForce_New(_particle, e_z, weightDiff_FaceZ);

  //Add force to cell faces
  m_cellFacesX[_cellIndex]->m_deviatoricForce+=force_FaceX;
  m_cellFacesY[_cellIndex]->m_deviatoricForce+=force_FaceY;
  m_cellFacesZ[_cellIndex]->m_deviatoricForce+=force_FaceZ;

  //Get mass of faces
  float mass_FaceX=m_cellFacesX[_cellIndex]->m_mass;
  float mass_FaceY=m_cellFacesY[_cellIndex]->m_mass;
  float mass_FaceZ=m_cellFacesZ[_cellIndex]->m_mass;

  //Calculate B Components
  float BComponentX=calcBComponent_DeviatoricVelocity_New(_particle, e_x, _weightX, force_FaceX, mass_FaceX);
  float BComponentY=calcBComponent_DeviatoricVelocity_New(_particle, e_y, _weightY, force_FaceY, mass_FaceY);
  float BComponentZ=calcBComponent_DeviatoricVelocity_New(_particle, e_z, _weightZ, force_FaceZ, mass_FaceZ);

  //Add contributions to B vectors
  m_Bvector_deviatoric_X(_cellIndex)+=BComponentX;
  m_Bvector_deviatoric_Y(_cellIndex)+=BComponentY;
  m_Bvector_deviatoric_Z(_cellIndex)+=BComponentZ;


  //Loop over cells the particle will contribute to again to get row components for A matrices

  //To calc position of particle, need origin of grid corner, not centre of first grid cell.
  float halfCellSize=m_cellSize/2.0;
  Eigen::Vector3f gridEdgePosition=m_origin;
  gridEdgePosition(0)-=halfCellSize;
  gridEdgePosition(1)-=halfCellSize;
  gridEdgePosition(2)-=halfCellSize;

  //Get particle position in grid
  Eigen::Vector3i particleIndex=MathFunctions::getParticleGridCell(particlePosition, m_cellSize, gridEdgePosition);

  //Loop over i+-2, j+-2, k+-2 to get cells that particle will contribute to
  int iParticle=particleIndex(0);
  int jParticle=particleIndex(1);
  int kParticle=particleIndex(2);

  for (int kIndex_neigh=(kParticle-2); kIndex_neigh<(kParticle+4); kIndex_neigh++)
  {
    for (int jIndex_neigh=(jParticle-2); jIndex_neigh<(jParticle+4); jIndex_neigh++)
    {
      for (int iIndex_neigh=(iParticle-2); iIndex_neigh<(iParticle+4); iIndex_neigh++)
      {
        //Check that not in outer cells or outside grid
        if (iIndex_neigh>=0 && iIndex_neigh<=(m_noCells-1) && jIndex_neigh>=0 && jIndex_neigh<=(m_noCells-1) && kIndex_neigh>=0 && kIndex_neigh<=(m_noCells-1))
        {
          //Get cell index
          int cellIndex_neigh=MathFunctions::getVectorIndex(iIndex_neigh, jIndex_neigh, kIndex_neigh, m_noCells);

          Eigen::Vector3f weightDiff_FaceX_neighbour;
          Eigen::Vector3f weightDiff_FaceY_neighbour;
          Eigen::Vector3f weightDiff_FaceZ_neighbour;

          //Get differentiated weights
          calcWeight_cubicBSpline_Diff(particlePosition, _iIndex, _jIndex, _kIndex, weightDiff_FaceX_neighbour, weightDiff_FaceY_neighbour, weightDiff_FaceZ_neighbour);

          //Calculate A components
          float AComponentX=calcAComponent_DeviatoricVelocity_New(_particle, weightDiff_FaceX, weightDiff_FaceX_neighbour, e_x);
          float AComponentY=calcAComponent_DeviatoricVelocity_New(_particle, weightDiff_FaceY, weightDiff_FaceY_neighbour, e_y);
          float AComponentZ=calcAComponent_DeviatoricVelocity_New(_particle, weightDiff_FaceZ, weightDiff_FaceZ_neighbour, e_z);

          m_Amatrix_deviatoric_X(_cellIndex, cellIndex_neigh)+=AComponentX;
          m_Amatrix_deviatoric_Y(_cellIndex, cellIndex_neigh)+=AComponentY;
          m_Amatrix_deviatoric_Z(_cellIndex, cellIndex_neigh)+=AComponentZ;

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


  //Addition so particle and grid is in same reference.
  //Origin of grid is at [-halfCell, -halfCell, -halfCell] in the particle reference system
  _particlePosition(0)+=halfCellSize;
  _particlePosition(1)+=halfCellSize;
  _particlePosition(2)+=halfCellSize;

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


//  //Centre
//  float dNx_cubicBS_Centre=MathFunctions::calcCubicBSpline_Diff(centrePosDiff(0)/m_cellSize);
//  float dNy_cubicBS_Centre=MathFunctions::calcCubicBSpline_Diff(centrePosDiff(1)/m_cellSize);
//  float dNz_cubicBS_Centre=MathFunctions::calcCubicBSpline_Diff(centrePosDiff(2)/m_cellSize);

//  o_weightCentre(0)=dNx_cubicBS_Centre*NyCentre_cubicBS*NzCentre_cubicBS;
//  o_weightCentre(1)=dNy_cubicBS_Centre*NxCentre_cubicBS*NzCentre_cubicBS;
//  o_weightCentre(2)=dNz_cubicBS_Centre*NxCentre_cubicBS*NyCentre_cubicBS;

//  o_weightCentre*=(1.0/m_cellSize); ///Not sure about this part?
////  o_weightCentre*=(-1.0/m_cellSize); ///Not sure about this part?

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

  o_weightFaceX(0)=dNx_cubicBS_FaceY*NyCentre_cubicBS*NzCentre_cubicBS;
  o_weightFaceX(1)=dNy_cubicBS_FaceY*NxCentre_cubicBS*NzCentre_cubicBS;
  o_weightFaceX(2)=dNz_cubicBS_FaceY*NxCentre_cubicBS*NyCentre_cubicBS;

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
