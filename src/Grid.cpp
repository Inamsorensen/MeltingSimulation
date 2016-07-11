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

  //Initialise time step to zero
  m_dt=0.0;

  //Initialise surrounding temperatures to zero
  m_ambientTemperature=0.0;
  m_heatSourceTemperature=0.0;


  m_cellCentres.reserve(pow(m_noCells,3));
  m_cellFacesX.reserve(pow(m_noCells,3));
  m_cellFacesY.reserve(pow(m_noCells,3));
  m_cellFacesZ.reserve(pow(m_noCells,3));

  //Setup cell lists
  for (int k=0; k<m_noCells; k++)
  {
    for (int j=0; j<m_noCells; j++)
    {
      for (int i=0; i<m_noCells; i++)
      {
        CellCentre* cellCentre=new CellCentre();
        CellFace* cellFaceX=new CellFace();
        CellFace* cellFaceY=new CellFace();
        CellFace* cellFaceZ=new CellFace();

        cellCentre->m_iIndex=i;
        cellCentre->m_jIndex=j;
        cellCentre->m_kIndex=k;

        cellFaceX->m_iIndex=i;
        cellFaceX->m_jIndex=j;
        cellFaceX->m_kIndex=k;

        cellFaceY->m_iIndex=i;
        cellFaceY->m_jIndex=j;
        cellFaceY->m_kIndex=k;

        cellFaceZ->m_iIndex=i;
        cellFaceZ->m_jIndex=j;
        cellFaceZ->m_kIndex=k;

        m_cellCentres.push_back(cellCentre);
        m_cellFacesX.push_back(cellFaceX);
        m_cellFacesY.push_back(cellFaceY);
        m_cellFacesZ.push_back(cellFaceZ);
      }
    }
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

void Grid::setSurroundingTemperatures(float _ambientTemp, float _heatSourceTemp)
{
  m_ambientTemperature=_ambientTemp+273.0;
  m_heatSourceTemperature=_heatSourceTemp+273.0;
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::update(float _dt, Emitter* _emitter, bool _isFirstStep)
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

  //Set m_dt=_dt
  m_dt=_dt;

  //Clear InterpolationData for each grid cell so all empty before start adding particles
  clearCellData();

  //findParticleInCell - need to find out which particles are in which cells and their respective interp weight
  findParticleInCell(_emitter);

  //Transfer particle data to grid
  transferParticleData(_emitter);

  //If first step calculate particle density during this loop as well
  if (_isFirstStep)
  {
    calcInitialParticleVolumes(_emitter);
  }

  //Classify cells
  classifyCells();


  //Calculate force

  //Calculate b

  //Set boundary velocity to be used in b
  setBoundaryVelocity();

  //Implicit integration to deviatoric velocity

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::clearCellData()
{
  /* Outline
  --------------------------------------------------------------------------------------------------------------
  Loop over all cells
    Clear m_interpolationData
    Set all variables to zero
    Set collision state of cell centre to colliding and faces to interior
  --------------------------------------------------------------------------------------------------------------
  */

  for (int cellIndex=0; cellIndex<pow(m_noCells,3); cellIndex++)
  {
    //Clear list of interpolation data including particle pointers.
    m_cellCentres[cellIndex]->m_interpolationData.clear();
    m_cellFacesX[cellIndex]->m_interpolationData.clear();
    m_cellFacesY[cellIndex]->m_interpolationData.clear();
    m_cellFacesZ[cellIndex]->m_interpolationData.clear();

    //Reset cell centre values to zero
    m_cellCentres[cellIndex]->m_mass=0.0;
    m_cellCentres[cellIndex]->m_detDeformationGrad=0.0;
    m_cellCentres[cellIndex]->m_detDeformationGradElastic=0.0;
    m_cellCentres[cellIndex]->m_detDeformationGradPlastic=0.0;
    m_cellCentres[cellIndex]->m_heatCapacity=0.0;
    m_cellCentres[cellIndex]->m_temperature=0.0;
    m_cellCentres[cellIndex]->m_lameLambdaInverse=0.0;
    m_cellCentres[cellIndex]->m_temperature=0.0;
    m_cellCentres[cellIndex]->m_state=State::Colliding;

    //Reset cell face X values to zero
    m_cellFacesX[cellIndex]->m_mass=0.0;
    m_cellFacesX[cellIndex]->m_deviatoricForce=0.0;
    m_cellFacesX[cellIndex]->m_velocity=0.0;
    m_cellFacesX[cellIndex]->m_heatConductivity=0.0;
    m_cellFacesX[cellIndex]->m_state=State::Interior;

    //Reset cell face Y values to zero
    m_cellFacesY[cellIndex]->m_mass=0.0;
    m_cellFacesY[cellIndex]->m_deviatoricForce=0.0;
    m_cellFacesY[cellIndex]->m_velocity=0.0;
    m_cellFacesY[cellIndex]->m_heatConductivity=0.0;
    m_cellFacesY[cellIndex]->m_state=State::Interior;

    //Reset cell face Z values to zero
    m_cellFacesZ[cellIndex]->m_mass=0.0;
    m_cellFacesZ[cellIndex]->m_deviatoricForce=0.0;
    m_cellFacesZ[cellIndex]->m_velocity=0.0;
    m_cellFacesZ[cellIndex]->m_heatConductivity=0.0;
    m_cellFacesZ[cellIndex]->m_state=State::Interior;

  }
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::findParticleInCell(Emitter* _emitter)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------
   Loop over all particles - This is where parallel should be inserted I think
   {
     Find position of particle in grid using getParticleGridCell
     This gives vector of i,j,k for cell

     Get neighbours i+-2, j+-2, k+-2. Ie loop over these - NB! watch out so don't check cells outside grid
     Actually needs to be i-1 and i+3 for faces, hence set -2 and +3 for now

     Pass in cell i,j,k to calcInterpolationWeights

    }
  ------------------------------------------------------------------------------------------------------
  */

  std::vector<Particle*>* particleListPtr=_emitter->getParticlesList();

  //To calc position of particle, need origin of grid ege, not centre of first grid cell, as this is how its
  //defined in Houdini/import file
  float halfCellSize=m_cellSize/2.0;
  Eigen::Vector3f gridEdgePosition=m_origin;
  gridEdgePosition(0)-=halfCellSize;
  gridEdgePosition(1)-=halfCellSize;
  gridEdgePosition(2)-=halfCellSize;

  for (int particleItr=0; particleItr<_emitter->getNoParticles(); particleItr++)
  {
    Eigen::Vector3f particlePosition=particleListPtr->at(particleItr)->getPosition();
    Eigen::Vector3i particleIndex=MathFunctions::getParticleGridCell(particlePosition, m_cellSize, gridEdgePosition);

    //Loop over i+-2, j+-2, k+-2
    int iParticle=particleIndex(0);
    int jParticle=particleIndex(1);
    int kParticle=particleIndex(2);

    for (int k=(kParticle-2); k<(kParticle+4); k++)
    {
      for (int j=(jParticle-2); j<(jParticle+4); j++)
      {
        for (int i=(iParticle-2); i<(iParticle+4); i++)
        {
          //Check that not in outer cells or outside grid
          if (i>=0 && i<=(m_noCells-1) && j>=0 && j<=(m_noCells-1) && k>=0 && k<=(m_noCells-1))
          {
            //Calculate interpolation weight and store particle if unlike zero
            calcInterpolationWeights(particleListPtr->at(particleItr), i, j, k);
          }
        }
      }
    }
  }
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

  //Cell position
  float xPos=(_i*m_cellSize)+m_origin(0);
  float yPos=(_j*m_cellSize)+m_origin(1);
  float zPos=(_k*m_cellSize)+m_origin(2);

  //Position vectors for centre and faces
  Eigen::Vector3f centreVector(xPos, yPos, zPos);

  float halfCellSize=m_cellSize/2.0;
  Eigen::Vector3f faceXVector(xPos-halfCellSize, yPos, zPos);
  Eigen::Vector3f faceYVector(xPos, yPos-halfCellSize, zPos);
  Eigen::Vector3f faceZVector(xPos, yPos, zPos-halfCellSize);

  Eigen::Vector3f particlePosition=_particle->getPosition();

  //Calculate posDifference for each face and cell centre
  Eigen::Vector3f centrePosDiff=particlePosition-centreVector;
  Eigen::Vector3f faceXPosDiff=particlePosition-faceXVector;
  Eigen::Vector3f faceYPosDiff=particlePosition-faceYVector;
  Eigen::Vector3f faceZPosDiff=particlePosition-faceZVector;

  //Need to calculate weights for each of these
  //Check whether worth calculating all?

  //Centre
  float NxCentre_cubicBS=MathFunctions::calcCubicBSpline(centrePosDiff(0)/m_cellSize);
  float NyCentre_cubicBS=MathFunctions::calcCubicBSpline(centrePosDiff(1)/m_cellSize);
  float NzCentre_cubicBS=MathFunctions::calcCubicBSpline(centrePosDiff(2)/m_cellSize);
  float NCentre_cubicBS=NxCentre_cubicBS*NyCentre_cubicBS*NzCentre_cubicBS;

  //FaceX
  float NxFaceX_cubicBS=MathFunctions::calcCubicBSpline(faceXPosDiff(0)/m_cellSize);
  float NyFaceX_cubicBS=MathFunctions::calcCubicBSpline(faceXPosDiff(1)/m_cellSize);
  float NzFaceX_cubicBS=MathFunctions::calcCubicBSpline(faceXPosDiff(2)/m_cellSize);
  float NFaceX_cubicBS=NxFaceX_cubicBS*NyFaceX_cubicBS*NzFaceX_cubicBS;

  //FaceY
  float NxFaceY_cubicBS=MathFunctions::calcCubicBSpline(faceYPosDiff(0)/m_cellSize);
  float NyFaceY_cubicBS=MathFunctions::calcCubicBSpline(faceYPosDiff(1)/m_cellSize);
  float NzFaceY_cubicBS=MathFunctions::calcCubicBSpline(faceYPosDiff(2)/m_cellSize);
  float NFaceY_cubicBS=NxFaceY_cubicBS*NyFaceY_cubicBS*NzFaceY_cubicBS;

  //FaceZ
  float NxFaceZ_cubicBS=MathFunctions::calcCubicBSpline(faceZPosDiff(0)/m_cellSize);
  float NyFaceZ_cubicBS=MathFunctions::calcCubicBSpline(faceZPosDiff(1)/m_cellSize);
  float NzFaceZ_cubicBS=MathFunctions::calcCubicBSpline(faceZPosDiff(2)/m_cellSize);
  float NFaceZ_cubicBS=NxFaceZ_cubicBS*NyFaceZ_cubicBS*NzFaceZ_cubicBS;

  //Check whether worth keep going, ie. if cubicBS are non-zero
  //NB! Might need to check if smaller than smallest value difference

  int cellListIndex=MathFunctions::getVectorIndex(_i, _j, _k, m_noCells);
  //Centre
  if (NCentre_cubicBS!=0)
  {
    //Create interpolation data pointer
    InterpolationData* newInterpolationData= new InterpolationData;

    //Store particle
    newInterpolationData->m_particle=_particle;

    //Store cubicBSpline
    newInterpolationData->m_cubicBSpline=NCentre_cubicBS;

    //Calculate and store cubicBSpline differentiated or nabla*weight for cubicBS weight
    Eigen::Vector3f cubicBS_Diff;

    float dNx_cubicBS=MathFunctions::calcCubicBSpline_Diff(centrePosDiff(0)/m_cellSize);
    float dNy_cubicBS=MathFunctions::calcCubicBSpline_Diff(centrePosDiff(1)/m_cellSize);
    float dNz_cubicBS=MathFunctions::calcCubicBSpline_Diff(centrePosDiff(2)/m_cellSize);

    cubicBS_Diff(0)=dNx_cubicBS*NyCentre_cubicBS*NzCentre_cubicBS;
    cubicBS_Diff(1)=dNy_cubicBS*NxCentre_cubicBS*NzCentre_cubicBS;
    cubicBS_Diff(2)=dNz_cubicBS*NxCentre_cubicBS*NyCentre_cubicBS;

    cubicBS_Diff*=(1.0/m_cellSize); ///Not sure about this part?
    newInterpolationData->m_cubicBSpline_Diff=cubicBS_Diff;


    //Calculate and store Tight Quadratic stencil
    float NxCentre_quadS=MathFunctions::calcTightQuadraticStencil(centrePosDiff(0)/m_cellSize);
    float NyCentre_quadS=MathFunctions::calcTightQuadraticStencil(centrePosDiff(1)/m_cellSize);
    float NzCentre_quadS=MathFunctions::calcTightQuadraticStencil(centrePosDiff(2)/m_cellSize);
    float NCentre_quadS=NxCentre_quadS*NyCentre_quadS*NzCentre_quadS;
    newInterpolationData->m_tightQuadStencil=NCentre_quadS;


    //Calculate and store Tight Quadratic stencil differentiated or nabla*weight for quadratic stencil
    Eigen::Vector3f quadS_Diff;

    float dNx_quadS=MathFunctions::calcTightQuadraticStencil_Diff(centrePosDiff(0)/m_cellSize);
    float dNy_quadS=MathFunctions::calcTightQuadraticStencil_Diff(centrePosDiff(1)/m_cellSize);
    float dNz_quadS=MathFunctions::calcTightQuadraticStencil_Diff(centrePosDiff(2)/m_cellSize);

    quadS_Diff(0)=dNx_quadS*NyCentre_quadS*NzCentre_quadS;
    quadS_Diff(1)=dNy_quadS*NxCentre_quadS*NzCentre_quadS;
    quadS_Diff(2)=dNz_quadS*NxCentre_quadS*NyCentre_quadS;

    quadS_Diff*=(1.0/m_cellSize); ///Not sure about this part?
    newInterpolationData->m_tightQuadStencil_Diff=quadS_Diff;


    //Store interpolation data
    m_cellCentres[cellListIndex]->m_interpolationData.push_back(newInterpolationData);


  }
  //FaceX
  if (NFaceX_cubicBS!=0)
  {
    //Create interpolation data pointer
    InterpolationData* newInterpolationData= new InterpolationData;

    //Store particle
    newInterpolationData->m_particle=_particle;

    //Store cubicBSpline
    newInterpolationData->m_cubicBSpline=NFaceX_cubicBS;

    //Calculate and store cubicBSpline differentiated or nabla*weight for cubicBS weight
    Eigen::Vector3f cubicBS_Diff;

    float dNx_cubicBS=MathFunctions::calcCubicBSpline_Diff(faceXPosDiff(0)/m_cellSize);
    float dNy_cubicBS=MathFunctions::calcCubicBSpline_Diff(faceXPosDiff(1)/m_cellSize);
    float dNz_cubicBS=MathFunctions::calcCubicBSpline_Diff(faceXPosDiff(2)/m_cellSize);

    cubicBS_Diff(0)=dNx_cubicBS*NyFaceX_cubicBS*NzFaceX_cubicBS;
    cubicBS_Diff(1)=dNy_cubicBS*NxFaceX_cubicBS*NzFaceX_cubicBS;
    cubicBS_Diff(2)=dNz_cubicBS*NxFaceX_cubicBS*NyFaceX_cubicBS;

    cubicBS_Diff*=(1.0/m_cellSize); ///Not sure about this part?
    newInterpolationData->m_cubicBSpline_Diff=cubicBS_Diff;


    //Calculate and store Tight Quadratic stencil
    float NxFaceX_quadS=MathFunctions::calcTightQuadraticStencil(faceXPosDiff(0)/m_cellSize);
    float NyFaceX_quadS=MathFunctions::calcTightQuadraticStencil(faceXPosDiff(1)/m_cellSize);
    float NzFaceX_quadS=MathFunctions::calcTightQuadraticStencil(faceXPosDiff(2)/m_cellSize);
    float NFaceX_quadS=NxFaceX_quadS*NyFaceX_quadS*NzFaceX_quadS;
    newInterpolationData->m_tightQuadStencil=NFaceX_quadS;


    //Calculate and store Tight Quadratic stencil differentiated or nabla*weight for quadratic stencil
    Eigen::Vector3f quadS_Diff;

    float dNx_quadS=MathFunctions::calcTightQuadraticStencil_Diff(faceXPosDiff(0)/m_cellSize);
    float dNy_quadS=MathFunctions::calcTightQuadraticStencil_Diff(faceXPosDiff(1)/m_cellSize);
    float dNz_quadS=MathFunctions::calcTightQuadraticStencil_Diff(faceXPosDiff(2)/m_cellSize);

    quadS_Diff(0)=dNx_quadS*NyFaceX_quadS*NzFaceX_quadS;
    quadS_Diff(1)=dNy_quadS*NxFaceX_quadS*NzFaceX_quadS;
    quadS_Diff(2)=dNz_quadS*NxFaceX_quadS*NyFaceX_quadS;

    quadS_Diff*=(1.0/m_cellSize); ///Not sure about this part?
    newInterpolationData->m_tightQuadStencil_Diff=quadS_Diff;


    //Store interpolation data
    m_cellFacesX[cellListIndex]->m_interpolationData.push_back(newInterpolationData);
  }
  //FaceY
  if (NFaceY_cubicBS!=0)
  {
    //Create interpolation data pointer
    InterpolationData* newInterpolationData= new InterpolationData;

    //Store particle
    newInterpolationData->m_particle=_particle;

    //Store cubicBSpline
    newInterpolationData->m_cubicBSpline=NFaceY_cubicBS;

    //Calculate and store cubicBSpline differentiated or nabla*weight for cubicBS weight
    Eigen::Vector3f cubicBS_Diff;

    float dNx_cubicBS=MathFunctions::calcCubicBSpline_Diff(faceYPosDiff(0)/m_cellSize);
    float dNy_cubicBS=MathFunctions::calcCubicBSpline_Diff(faceYPosDiff(1)/m_cellSize);
    float dNz_cubicBS=MathFunctions::calcCubicBSpline_Diff(faceYPosDiff(2)/m_cellSize);

    cubicBS_Diff(0)=dNx_cubicBS*NyFaceY_cubicBS*NzFaceY_cubicBS;
    cubicBS_Diff(1)=dNy_cubicBS*NxFaceY_cubicBS*NzFaceY_cubicBS;
    cubicBS_Diff(2)=dNz_cubicBS*NxFaceY_cubicBS*NyFaceY_cubicBS;

    cubicBS_Diff*=(1.0/m_cellSize); ///Not sure about this part?
    newInterpolationData->m_cubicBSpline_Diff=cubicBS_Diff;


    //Calculate and store Tight Quadratic stencil
    float NxFaceY_quadS=MathFunctions::calcTightQuadraticStencil(faceYPosDiff(0)/m_cellSize);
    float NyFaceY_quadS=MathFunctions::calcTightQuadraticStencil(faceYPosDiff(1)/m_cellSize);
    float NzFaceY_quadS=MathFunctions::calcTightQuadraticStencil(faceYPosDiff(2)/m_cellSize);
    float NFaceY_quadS=NxFaceY_quadS*NyFaceY_quadS*NzFaceY_quadS;
    newInterpolationData->m_tightQuadStencil=NFaceY_quadS;


    //Calculate and store Tight Quadratic stencil differentiated or nabla*weight for quadratic stencil
    Eigen::Vector3f quadS_Diff;

    float dNx_quadS=MathFunctions::calcTightQuadraticStencil_Diff(faceYPosDiff(0)/m_cellSize);
    float dNy_quadS=MathFunctions::calcTightQuadraticStencil_Diff(faceYPosDiff(1)/m_cellSize);
    float dNz_quadS=MathFunctions::calcTightQuadraticStencil_Diff(faceYPosDiff(2)/m_cellSize);

    quadS_Diff(0)=dNx_quadS*NyFaceY_quadS*NzFaceY_quadS;
    quadS_Diff(1)=dNy_quadS*NxFaceY_quadS*NzFaceY_quadS;
    quadS_Diff(2)=dNz_quadS*NxFaceY_quadS*NyFaceY_quadS;

    quadS_Diff*=(1.0/m_cellSize); ///Not sure about this part?
    newInterpolationData->m_tightQuadStencil_Diff=quadS_Diff;


    //Store interpolation data
    m_cellFacesY[cellListIndex]->m_interpolationData.push_back(newInterpolationData);
  }
  //FaceZ
  if (NFaceZ_cubicBS!=0)
  {
    //Create interpolation data pointer
    InterpolationData* newInterpolationData= new InterpolationData;

    //Store particle
    newInterpolationData->m_particle=_particle;

    //Store cubicBSpline
    newInterpolationData->m_cubicBSpline=NFaceZ_cubicBS;

    //Calculate and store cubicBSpline differentiated or nabla*weight for cubicBS weight
    Eigen::Vector3f cubicBS_Diff;

    float dNx_cubicBS=MathFunctions::calcCubicBSpline_Diff(faceZPosDiff(0)/m_cellSize);
    float dNy_cubicBS=MathFunctions::calcCubicBSpline_Diff(faceZPosDiff(1)/m_cellSize);
    float dNz_cubicBS=MathFunctions::calcCubicBSpline_Diff(faceZPosDiff(2)/m_cellSize);

    cubicBS_Diff(0)=dNx_cubicBS*NyFaceZ_cubicBS*NzFaceZ_cubicBS;
    cubicBS_Diff(1)=dNy_cubicBS*NxFaceZ_cubicBS*NzFaceZ_cubicBS;
    cubicBS_Diff(2)=dNz_cubicBS*NxFaceZ_cubicBS*NyFaceZ_cubicBS;

    cubicBS_Diff*=(1.0/m_cellSize); ///Not sure about this part?
    newInterpolationData->m_cubicBSpline_Diff=cubicBS_Diff;


    //Calculate and store Tight Quadratic stencil
    float NxFaceZ_quadS=MathFunctions::calcTightQuadraticStencil(faceZPosDiff(0)/m_cellSize);
    float NyFaceZ_quadS=MathFunctions::calcTightQuadraticStencil(faceZPosDiff(1)/m_cellSize);
    float NzFaceZ_quadS=MathFunctions::calcTightQuadraticStencil(faceZPosDiff(2)/m_cellSize);
    float NFaceZ_quadS=NxFaceZ_quadS*NyFaceZ_quadS*NzFaceZ_quadS;
    newInterpolationData->m_tightQuadStencil=NFaceZ_quadS;


    //Calculate and store Tight Quadratic stencil differentiated or nabla*weight for quadratic stencil
    Eigen::Vector3f quadS_Diff;

    float dNx_quadS=MathFunctions::calcTightQuadraticStencil_Diff(faceZPosDiff(0)/m_cellSize);
    float dNy_quadS=MathFunctions::calcTightQuadraticStencil_Diff(faceZPosDiff(1)/m_cellSize);
    float dNz_quadS=MathFunctions::calcTightQuadraticStencil_Diff(faceZPosDiff(2)/m_cellSize);

    quadS_Diff(0)=dNx_quadS*NyFaceZ_quadS*NzFaceZ_quadS;
    quadS_Diff(1)=dNy_quadS*NxFaceZ_quadS*NzFaceZ_quadS;
    quadS_Diff(2)=dNz_quadS*NxFaceZ_quadS*NyFaceZ_quadS;

    quadS_Diff*=(1.0/m_cellSize); ///Not sure about this part?
    newInterpolationData->m_tightQuadStencil_Diff=quadS_Diff;


    //Store interpolation data
    m_cellFacesZ[cellListIndex]->m_interpolationData.push_back(newInterpolationData);
  }




}

//----------------------------------------------------------------------------------------------------------------------

void Grid::transferParticleData(Emitter* _emitter)
{
  /* Outline
  -----------------------------------------------------------------------------------------------------
  Loop done in update - could be moved to here and remove cellIndex input

    Check that m_InterpolationData is not empty - If it is, set everything to zero?

    For each particle in the list calculate all variables, i for cell faces i={x,y,z} and c for cell centre
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

  -----------------------------------------------------------------------------------------------------
  */

  for (int cellIndex=0; cellIndex<pow(m_noCells, 3); cellIndex++)
  {
    //Face X
    //Check that non-empty, ie. that it has particles in it
    int noParticles_CellFaceX=m_cellFacesX[cellIndex]->m_interpolationData.size();
    if (noParticles_CellFaceX!=0)
    {
      for (int particleIterator=0; particleIterator<noParticles_CellFaceX; particleIterator++)
      {
        //Get interpolation weight
        float weight=m_cellFacesX[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

        //Get particle data
        float mass=0.0;
        Eigen::Vector3f velocity;
        Phase phase=Phase::Solid;
        m_cellFacesX[cellIndex]->m_interpolationData[particleIterator]->m_particle->getParticleData_CellFace(mass, velocity, phase);
        float velocityX=velocity(0);

        //Add to cell face data
        m_cellFacesX[cellIndex]->m_mass+=(weight*mass);
        m_cellFacesX[cellIndex]->m_velocity+=((weight*mass)*velocityX);

        //Find heat conductivity depending on phase
        float heatConductivity=0.0;
        if (phase==Phase::Solid)
        {
          heatConductivity=_emitter->m_heatConductivitySolid;
        }
        else
        {
          heatConductivity=_emitter->m_heatConductivityFluid;
        }
        m_cellFacesX[cellIndex]->m_heatConductivity+=((weight*mass)*heatConductivity);

      }

      //Multiply data by 1/m_{i}
      m_cellFacesX[cellIndex]->m_velocity*=(1.0/m_cellFacesX[cellIndex]->m_mass);
      m_cellFacesX[cellIndex]->m_heatConductivity*=(1.0/m_cellFacesX[cellIndex]->m_mass);

    }

    //Face Y
    int noParticles_CellFaceY=m_cellFacesY[cellIndex]->m_interpolationData.size();
    if (noParticles_CellFaceY!=0)
    {
      for (int particleIterator=0; particleIterator<noParticles_CellFaceY; particleIterator++)
      {
        //Get interpolation weight
        float weight=m_cellFacesY[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

        //Get particle data
        float mass=0.0;
        Eigen::Vector3f velocity;
        Phase phase=Phase::Solid;
        m_cellFacesY[cellIndex]->m_interpolationData[particleIterator]->m_particle->getParticleData_CellFace(mass, velocity, phase);
        float velocityY=velocity(1);

        //Add to cell face data
        m_cellFacesY[cellIndex]->m_mass+=(weight*mass);
        m_cellFacesY[cellIndex]->m_velocity+=((weight*mass)*velocityY);

        //Find heat conductivity depending on phase
        float heatConductivity=0.0;
        if (phase==Phase::Solid)
        {
          heatConductivity=_emitter->m_heatConductivitySolid;
        }
        else
        {
          heatConductivity=_emitter->m_heatConductivityFluid;
        }
        m_cellFacesY[cellIndex]->m_heatConductivity+=((weight*mass)*heatConductivity);

      }

      //Multiply data by 1/m_{i}
      m_cellFacesY[cellIndex]->m_velocity*=(1.0/m_cellFacesY[cellIndex]->m_mass);
      m_cellFacesY[cellIndex]->m_heatConductivity*=(1.0/m_cellFacesY[cellIndex]->m_mass);

    }

    //Face Z
    int noParticles_CellFaceZ=m_cellFacesZ[cellIndex]->m_interpolationData.size();
    if (noParticles_CellFaceZ!=0)
    {
      for (int particleIterator=0; particleIterator<noParticles_CellFaceZ; particleIterator++)
      {
        //Get interpolation weight
        float weight=m_cellFacesZ[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

        //Get particle data
        float mass=0.0;
        Eigen::Vector3f velocity;
        Phase phase=Phase::Solid;
        m_cellFacesZ[cellIndex]->m_interpolationData[particleIterator]->m_particle->getParticleData_CellFace(mass, velocity, phase);
        float velocityZ=velocity(2);

        //Add to cell face data
        m_cellFacesZ[cellIndex]->m_mass+=(weight*mass);
        m_cellFacesZ[cellIndex]->m_velocity+=((weight*mass)*velocityZ);

        //Find heat conductivity depending on phase
        float heatConductivity=0.0;
        if (phase==Phase::Solid)
        {
          heatConductivity=_emitter->m_heatConductivitySolid;
        }
        else
        {
          heatConductivity=_emitter->m_heatConductivityFluid;
        }
        m_cellFacesZ[cellIndex]->m_heatConductivity+=((weight*mass)*heatConductivity);

      }

      //Multiply data by 1/m_{i}
      m_cellFacesZ[cellIndex]->m_velocity*=(1.0/m_cellFacesZ[cellIndex]->m_mass);
      m_cellFacesZ[cellIndex]->m_heatConductivity*=(1.0/m_cellFacesZ[cellIndex]->m_mass);

    }

    //Cell centre
    int noParticles_CellCentre=m_cellCentres[cellIndex]->m_interpolationData.size();
    if (noParticles_CellCentre!=0)
    {
      for (int particleIterator=0; particleIterator<noParticles_CellCentre; particleIterator++)
      {
        //Get interpolation weight
        float weight=m_cellCentres[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

        //Get particle data
        float mass=0.0;
        float detDeformGrad=0.0;
        float detDeformGradElast=0.0;
        Phase phase=Phase::Solid;
        float temperature=0.0;
        float lameLambdaInverse=0.0;
        m_cellCentres[cellIndex]->m_interpolationData[particleIterator]->m_particle->getParticleData_CellCentre(mass, detDeformGrad, detDeformGradElast, phase, temperature, lameLambdaInverse);

        //Add to cell centre data
        m_cellCentres[cellIndex]->m_mass+=(weight*mass);
        m_cellCentres[cellIndex]->m_detDeformationGrad+=((weight*mass)*detDeformGrad);
        m_cellCentres[cellIndex]->m_detDeformationGradElastic+=((weight*mass)*detDeformGradElast);
        m_cellCentres[cellIndex]->m_temperature+=((weight*mass)*temperature);
        m_cellCentres[cellIndex]->m_lameLambdaInverse+=((weight*mass)*lameLambdaInverse);

        //Get heat capacity depending on phase
        float heatCapacity=0.0;
        if (phase==Phase::Solid)
        {
          heatCapacity=_emitter->m_heatCapacitySolid;
        }
        else
        {
          heatCapacity=_emitter->m_heatCapacityFluid;
        }
        m_cellCentres[cellIndex]->m_heatCapacity+=((weight*mass)*heatCapacity);

      }

      //Multiply data by 1/m_{c}
      m_cellCentres[cellIndex]->m_detDeformationGrad*=(1.0/m_cellCentres[cellIndex]->m_mass);
      m_cellCentres[cellIndex]->m_detDeformationGradElastic*=(1.0/m_cellCentres[cellIndex]->m_mass);
      m_cellCentres[cellIndex]->m_heatCapacity*=(1.0/m_cellCentres[cellIndex]->m_mass);
      m_cellCentres[cellIndex]->m_temperature*=(1.0/m_cellCentres[cellIndex]->m_mass);
      m_cellCentres[cellIndex]->m_lameLambdaInverse*=(1.0/m_cellCentres[cellIndex]->m_mass);

      //Calculate detDeformationGrad_Plastic, ie. J_{Pc}=J_{c}/J_{Ec}
      m_cellCentres[cellIndex]->m_detDeformationGradPlastic=m_cellCentres[cellIndex]->m_detDeformationGrad;
      m_cellCentres[cellIndex]->m_detDeformationGradPlastic*=(1.0/m_cellCentres[cellIndex]->m_detDeformationGradElastic);
    }
  }
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcInitialParticleVolumes(Emitter *_emitter)
{
  for (int cellIndex=0; cellIndex<pow(m_noCells, 3); cellIndex++)
  {
    //Cell volume
    float cellVolume=pow(m_cellSize,3);

    //Add grid cells contribution to particle density
    int noParticles_CellCentre=m_cellCentres[cellIndex]->m_interpolationData.size();

    //Get cell centre mass
    float mass=m_cellCentres[cellIndex]->m_mass;

    for (int particleIterator=0; particleIterator<noParticles_CellCentre; particleIterator++)
    {
      //Get cubicBSpline weight for cell centre i and particle particleIterator
      float weight=m_cellCentres[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

      //Add density from this cell to particle
      float density=(weight*mass)/cellVolume;
      m_cellCentres[cellIndex]->m_interpolationData[particleIterator]->m_particle->addParticleDensity(density);
//      std::cout<<"test\n";
    }
  }

  //Calculate particle volume
  int noParticles=_emitter->getNoParticles();
  std::vector<Particle*>* particleListPtr=_emitter->getParticlesList();
  for (int particleIterator=0; particleIterator<noParticles; particleIterator++)
  {
    particleListPtr->at(particleIterator)->calcInitialVolume();
  }
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
  for (int cellIndex=0; cellIndex<(pow(m_noCells, 3)); cellIndex++)
  {
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

  //Loop over all cells again to check which cell centres are collding
  //Seems inefficient.
  for (int cellIndex=0; cellIndex<(pow(m_noCells, 3)); cellIndex++)
  {
    //Find cell index.
    int iIndex=m_cellCentres[cellIndex]->m_iIndex;
    int jIndex=m_cellCentres[cellIndex]->m_jIndex;
    int kIndex=m_cellCentres[cellIndex]->m_kIndex;


    //Face X
    //Check if lower x face colliding
    if (m_cellFacesX[cellIndex]->m_state!=State::Colliding)
    {
      //Check whether empty or not
      int noParticlesInCell=m_cellCentres[cellIndex]->m_interpolationData.size();

      if (noParticlesInCell!=0)
      {
        m_cellCentres[cellIndex]->m_state=State::Interior;
      }
      else
      {
        m_cellCentres[cellIndex]->m_state=State::Empty;
        m_cellCentres[cellIndex]->m_temperature=m_ambientTemperature;
      }

      //Go to next cellIndex
      continue;
    }

    //Check upper x face as well
    if (iIndex<(m_noCells-1))
    {
      //Get index of cell with X face next to current cell
      int neighbourFaceIndex=MathFunctions::getVectorIndex(iIndex+1, jIndex, kIndex, m_noCells);

      if (m_cellFacesX[neighbourFaceIndex]->m_state!=State::Colliding)
      {
        //Check whether empty or not
        int noParticlesInCell=m_cellCentres[cellIndex]->m_interpolationData.size();

        if (noParticlesInCell!=0)
        {
          m_cellCentres[cellIndex]->m_state=State::Interior;
        }
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
      int noParticlesInCell=m_cellCentres[cellIndex]->m_interpolationData.size();

      if (noParticlesInCell!=0)
      {
        m_cellCentres[cellIndex]->m_state=State::Interior;
      }
      else
      {
        m_cellCentres[cellIndex]->m_state=State::Empty;
        m_cellCentres[cellIndex]->m_temperature=m_ambientTemperature;
      }

      //Go to next cellIndex
      continue;
    }

    //Check upper x face as well
    if (jIndex<(m_noCells-1))
    {
      //Get index of cell with Y face next to current cell
      int neighbourFaceIndex=MathFunctions::getVectorIndex(iIndex, jIndex+1, kIndex, m_noCells);

      if (m_cellFacesY[neighbourFaceIndex]->m_state!=State::Colliding)
      {
        //Check whether empty or not
        int noParticlesInCell=m_cellCentres[cellIndex]->m_interpolationData.size();

        if (noParticlesInCell!=0)
        {
          m_cellCentres[cellIndex]->m_state=State::Interior;
        }
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
      int noParticlesInCell=m_cellCentres[cellIndex]->m_interpolationData.size();

      if (noParticlesInCell!=0)
      {
        m_cellCentres[cellIndex]->m_state=State::Interior;
      }
      else
      {
        m_cellCentres[cellIndex]->m_state=State::Empty;
        m_cellCentres[cellIndex]->m_temperature=m_ambientTemperature;
      }

      //Go to next cellIndex
      continue;
    }

    //Check upper x face as well
    if (kIndex<(m_noCells-1))
    {
      //Get index of cell with Z face next to current cell
      int neighbourFaceIndex=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex+1, m_noCells);

      if (m_cellFacesZ[neighbourFaceIndex]->m_state!=State::Colliding)
      {
        //Check whether empty or not
        int noParticlesInCell=m_cellCentres[cellIndex]->m_interpolationData.size();

        if (noParticlesInCell!=0)
        {
          m_cellCentres[cellIndex]->m_state=State::Interior;
        }
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
    else if (m_cellCentres[cellIndex]->m_interpolationData.size()==0)
    {
      m_cellCentres[cellIndex]->m_temperature=m_ambientTemperature;
    }

  }
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::setBoundaryVelocity()
{
  /* Outline
  ---------------------------------------------------------------------------------------------------------------------
  Loop over all cells
    Check if faces are actually boundaries to colliding cell
    Done by checking whether current cell is colliding. If yes, then all faces are collision boundaries
    If current cell not collidining but faces are, check cells before them
      If cell before doesn't exist cause current cell is ijk=0, then set faces as collision boundaries
      Else check cell before in either face direction and if that cell is colliding, set face to collision boundary

    When face found to be collision boundary, set its velocity to zero for stick.
    This will need to change if object it's colliding with has a velocity.

  ---------------------------------------------------------------------------------------------------------------------
  */

  for (int cellIndex=0; cellIndex<(pow(m_noCells,3)); cellIndex++)
  {
    //Check whether current cell is colliding, if so set all faces to colliding
    if (m_cellCentres[cellIndex]->m_state==State::Colliding)
    {
      m_cellFacesX[cellIndex]->m_velocity=0.0;
      m_cellFacesY[cellIndex]->m_velocity=0.0;
      m_cellFacesZ[cellIndex]->m_velocity=0.0;
    }
    else
    {
      //Get cell index in ijk values
      int iIndex=m_cellCentres[cellIndex]->m_iIndex;
      int jIndex=m_cellCentres[cellIndex]->m_jIndex;
      int kIndex=m_cellCentres[cellIndex]->m_kIndex;

      //If current cell isn't colliding, then still need to check the cell before in each direction, unless i||j||k=0
      //FaceX
      if (m_cellFacesX[cellIndex]->m_state==State::Colliding)
      {
        //Set velocity to zero if iIndex=0
        if (iIndex==0)
        {
          m_cellFacesX[cellIndex]->m_velocity=0.0;
        }
        //Else check the cell before it in i direction
        else
        {
          //Get index of cell before in i direction
          int neighbourIndex=MathFunctions::getVectorIndex(iIndex-1, jIndex, kIndex, m_noCells);

          //Check if colliding
          if (m_cellCentres[neighbourIndex]->m_state==State::Colliding)
          {
            m_cellFacesX[cellIndex]->m_velocity=0.0;
          }
        }
      }

      //FaceY
      if (m_cellFacesY[cellIndex]->m_state==State::Colliding)
      {
        //Set velocity to zero if jIndex=0
        if (jIndex==0)
        {
          m_cellFacesY[cellIndex]->m_velocity=0.0;
        }
        //Else check the cell before it in j direction
        else
        {
          //Get index of cell before in j direction
          int neighbourIndex=MathFunctions::getVectorIndex(iIndex, jIndex-1, kIndex, m_noCells);

          //Check if colliding
          if (m_cellCentres[neighbourIndex]->m_state==State::Colliding)
          {
            m_cellFacesY[cellIndex]->m_velocity=0.0;
          }
        }
      }
      //FaceZ
      if (m_cellFacesZ[cellIndex]->m_state==State::Colliding)
      {
        //Set velocity to zero if kIndex=0
        if (kIndex==0)
        {
          m_cellFacesZ[cellIndex]->m_velocity=0.0;
        }
        //Else check the cell before it in k direction
        else
        {
          //Get index of cell before in k direction
          int neighbourIndex=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex-1, m_noCells);

          //Check if colliding
          if (m_cellCentres[neighbourIndex]->m_state==State::Colliding)
          {
            m_cellFacesZ[cellIndex]->m_velocity=0.0;
          }
        }
      }
    }
  }

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::findNoParticlesInCells(Emitter *_emitter, std::vector<int> &o_listParticleNo)
{
  /* Outline
  ---------------------------------------------------------------------------------------------------------------------
  Set up vector so large enough to contain all cells
  Initialise all members to zero

  Loop over all particles
    Find index of particle position
    Add to that cell in the vector

  ---------------------------------------------------------------------------------------------------------------------
  */

  //Initialise list
  //o_listParticleNo=new std::vector<int>(pow(m_noCells,3),0);

  //Get particle list
  std::vector<Particle*>* particleList=_emitter->getParticlesList();
  int noParticles=_emitter->getNoParticles();

  //Calculate position of grid edge as this is origin for particle positions
  float halfCellSize=m_cellSize/2.0;
  Eigen::Vector3f gridEdgePosition=m_origin;
  gridEdgePosition(0)-=halfCellSize;
  gridEdgePosition(1)-=halfCellSize;
  gridEdgePosition(2)-=halfCellSize;

  for (int i=0; i<noParticles; i++)
  {
    //Get grid cell index from particle position
    Eigen::Vector3f particlePosition=particleList->at(i)->getPosition();
    Eigen::Vector3i particleIndex=MathFunctions::getParticleGridCell(particlePosition, m_cellSize, gridEdgePosition);

    //Get vector index of the cell the particle is in
    int cellIndex=MathFunctions::getVectorIndex(particleIndex(0), particleIndex(1), particleIndex(2), m_noCells);

    //Increase the particle count for that cell
    o_listParticleNo.at(cellIndex)+=1;
  }

}

//----------------------------------------------------------------------------------------------------------------------

//void Grid::TEST_findParticleInCell(Emitter *_emitter)
//{
//  std::vector<Particle*>* particleListPtr=_emitter->getParticlesList();

//  for (int particleItr=0; particleItr<_emitter->getNoParticles(); particleItr++)
//  {
//    Eigen::Vector3f particlePosition=particleListPtr->at(particleItr)->getPosition();
//    //Eigen::Vector3i particleIndex=MathFunctions::getParticleGridCell(particlePosition, m_cellSize, m_origin);

//    //Loop over all grid cells
//    for (int k=0; k<m_noCells; k++)
//    {
//      for (int j=0; j<m_noCells; j++)
//      {
//        for (int i=0; i<m_noCells; i++)
//        {
//          //Cell position
//          float xPos=(i*m_cellSize)+m_origin(0);
//          float yPos=(j*m_cellSize)+m_origin(1);
//          float zPos=(k*m_cellSize)+m_origin(2);

//          //Position vectors for centre and faces
//          Eigen::Vector3f centreVector(xPos, yPos, zPos);

//          float halfCellSize=m_cellSize/2.0;
//          Eigen::Vector3f faceXVector(xPos-halfCellSize, yPos, zPos);
//          Eigen::Vector3f faceYVector(xPos, yPos-halfCellSize, zPos);
//          Eigen::Vector3f faceZVector(xPos, yPos, zPos-halfCellSize);

//          //Eigen::Vector3f particlePosition=_particle->getPosition();

//          //Calculate posDifference for each face and cell centre
//          Eigen::Vector3f centrePosDiff=particlePosition-centreVector;
//          Eigen::Vector3f faceXPosDiff=particlePosition-faceXVector;
//          Eigen::Vector3f faceYPosDiff=particlePosition-faceYVector;
//          Eigen::Vector3f faceZPosDiff=particlePosition-faceZVector;

//          //Need to calculate weights for each of these
//          //Check whether worth calculating all?

//          //Centre
//          float NxCentre_cubicBS=MathFunctions::calcCubicBSpline(centrePosDiff(0)/m_cellSize);
//          float NyCentre_cubicBS=MathFunctions::calcCubicBSpline(centrePosDiff(1)/m_cellSize);
//          float NzCentre_cubicBS=MathFunctions::calcCubicBSpline(centrePosDiff(2)/m_cellSize);
//          float NCentre_cubicBS=NxCentre_cubicBS*NyCentre_cubicBS*NzCentre_cubicBS;

//          //FaceX
//          float NxFaceX_cubicBS=MathFunctions::calcCubicBSpline(faceXPosDiff(0)/m_cellSize);
//          float NyFaceX_cubicBS=MathFunctions::calcCubicBSpline(faceXPosDiff(1)/m_cellSize);
//          float NzFaceX_cubicBS=MathFunctions::calcCubicBSpline(faceXPosDiff(2)/m_cellSize);
//          float NFaceX_cubicBS=NxFaceX_cubicBS*NyFaceX_cubicBS*NzFaceX_cubicBS;

//          //FaceY
//          float NxFaceY_cubicBS=MathFunctions::calcCubicBSpline(faceYPosDiff(0)/m_cellSize);
//          float NyFaceY_cubicBS=MathFunctions::calcCubicBSpline(faceYPosDiff(1)/m_cellSize);
//          float NzFaceY_cubicBS=MathFunctions::calcCubicBSpline(faceYPosDiff(2)/m_cellSize);
//          float NFaceY_cubicBS=NxFaceY_cubicBS*NyFaceY_cubicBS*NzFaceY_cubicBS;

//          //FaceZ
//          float NxFaceZ_cubicBS=MathFunctions::calcCubicBSpline(faceZPosDiff(0)/m_cellSize);
//          float NyFaceZ_cubicBS=MathFunctions::calcCubicBSpline(faceZPosDiff(1)/m_cellSize);
//          float NzFaceZ_cubicBS=MathFunctions::calcCubicBSpline(faceZPosDiff(2)/m_cellSize);
//          float NFaceZ_cubicBS=NxFaceZ_cubicBS*NyFaceZ_cubicBS*NzFaceZ_cubicBS;

//          //Check whether worth keep going, ie. if cubicBS are non-zero
//          //NB! Might need to check if smaller than smallest value difference

//          int cellListIndex=MathFunctions::getVectorIndex(i, j, k, m_noCells);

//          if (NCentre_cubicBS!=0)
//          {
//            InterpolationData* newInterpData=new InterpolationData;
//            newInterpData->m_particle=particleListPtr->at(particleItr);
//            newInterpData->m_cubicBSpline=NCentre_cubicBS;
//            m_cellCentres[cellListIndex]->m_interpolationData.push_back(newInterpData);
//          }

//          if (NFaceX_cubicBS!=0)
//          {
//            InterpolationData* newInterpData=new InterpolationData;
//            newInterpData->m_particle=particleListPtr->at(particleItr);
//            newInterpData->m_cubicBSpline=NFaceX_cubicBS;
//            m_cellFacesX[cellListIndex]->m_interpolationData.push_back(newInterpData);
//          }

//          if (NFaceY_cubicBS!=0)
//          {
//            InterpolationData* newInterpData=new InterpolationData;
//            newInterpData->m_particle=particleListPtr->at(particleItr);
//            newInterpData->m_cubicBSpline=NFaceY_cubicBS;
//            m_cellFacesY[cellListIndex]->m_interpolationData.push_back(newInterpData);
//          }

//          if (NFaceZ_cubicBS!=0)
//          {
//            InterpolationData* newInterpData=new InterpolationData;
//            newInterpData->m_particle=particleListPtr->at(particleItr);
//            newInterpData->m_cubicBSpline=NFaceZ_cubicBS;
//            m_cellFacesZ[cellListIndex]->m_interpolationData.push_back(newInterpData);
//          }
//        }
//      }
//    }
//  }
//}
