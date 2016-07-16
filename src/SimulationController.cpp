#include  <iostream>

#include "SimulationController.h"

//----------------------------------------------------------------------------------------------------------------------

SimulationController* SimulationController::m_instance=nullptr;

//----------------------------------------------------------------------------------------------------------------------

SimulationController::SimulationController()
{
  /// @brief Set up simulation parameters and create grid and emitter pointers

  //Initialise render parameters to zero. Set by separate function
  m_camera=nullptr;
  m_shaderName="";
  m_particleRadius=0.0;

  //Setup if don't read from file
  //Timer setup
  m_simTimeStep=0.01;
  m_elapsedTimeAfterFrame=0.0;
  m_noFrames=0;

  //Initialise total number of frames to zero
  m_totalNoFrames=0;

  //Grid setup
  Eigen::Vector3f boundingBoxPos;
  boundingBoxPos(0)=-0.5;
  boundingBoxPos(1)=-0.5;
  boundingBoxPos(2)=-0.5;
  m_boundingBoxPosition=boundingBoxPos;
  m_boundingBoxSize=1.0;
  m_noCells=16;

  //Particle setup
  m_noParticles=10000;
  m_particleMass=0.1;

  //Material setup
  m_lameMuConstant=1.0;
  m_lameLambdaConstant=1.0;
  m_hardnessCoefficient=1.0;
  m_compressionLimit=1.0;
  m_stretchLimit=1.0;

  m_heatCapacitySolid=1.0;
  m_heatCapacityFluid=1.0;
  m_heatConductivitySolid=1.0;
  m_heatConductivityFluid=1.0;
  m_latentHeat=1.0;
  m_freezingTemperature=0.0;

  //Set PIC FLIP contribution constants
  m_velocityContributionAlpha=0.95;
  m_temperatureContributionBeta=0.95;


  //Read in simulation parameters
  m_readFileName="../HoudiniFiles/particles.geo";
  readSimulationParameters();

  //Create emitter and particles
  m_emitter=new Emitter();
  m_emitter->setStrainConstants(m_lameMuConstant, m_lameLambdaConstant, m_compressionLimit, m_stretchLimit, m_hardnessCoefficient);
  m_emitter->setTemperatureConstants(m_heatCapacitySolid, m_heatCapacityFluid, m_heatConductivitySolid, m_heatConductivityFluid, m_latentHeat, m_freezingTemperature);
  setupParticles();

  //Create grid
  m_grid=Grid::createGrid(m_boundingBoxPosition, m_boundingBoxSize, m_noCells);
  m_grid->setSurroundingTemperatures(m_ambientTemperature, m_heatSourceTemperature);

  //Set grid as collision object for emitter
  float xMin=m_boundingBoxPosition(0);
  float yMin=m_boundingBoxPosition(1);
  float zMin=m_boundingBoxPosition(2);
  float xMax=xMin+m_boundingBoxSize;
  float yMax=yMin+m_boundingBoxSize;
  float zMax=zMin+m_boundingBoxSize;
  m_emitter->setCollisionObject(xMin, xMax, yMin, yMax, zMin, zMax);

  //Test min no particle in non-empty cells
  std::vector<int> listParticleNoInCells(pow(m_noCells,3),0);
  m_grid->findNoParticlesInCells(m_emitter, listParticleNoInCells);
  int minNoParticles=MathFunctions::findMinVectorValue(listParticleNoInCells);
  std::cout<<"The smallest number of particles in a non-empty cell is: "<<minNoParticles<<"\n";

}

//----------------------------------------------------------------------------------------------------------------------

SimulationController::~SimulationController()
{
  /// @brief Delete all pointers

  //Delete emitter and grid pointers
  delete m_emitter;
  delete m_grid;

  std::cout<<"Removing simulation controller\n";

}

//----------------------------------------------------------------------------------------------------------------------

SimulationController* SimulationController::instance()
{
  /// @brief Create simulation controller if doesn't exist, then return instance pointer

  if (m_instance==nullptr)
  {
    m_instance=new SimulationController();
  }

  return m_instance;
}

//----------------------------------------------------------------------------------------------------------------------

void SimulationController::setRenderParameters(ngl::Camera *_camera, std::string _shaderName)
{
  /// @brief Sets paramteres for rendering and the particle size
  /// @todo: Not sure what use the camera is. Will change as window is resized

  //Set camera. Not sure what use this is?
  m_camera=_camera;

  m_shaderName=_shaderName;

  //Particle size
  m_particleRadius=0.2;

  m_emitter->setRenderParameters(m_shaderName, m_particleRadius);

}

//----------------------------------------------------------------------------------------------------------------------

void SimulationController::readSimulationParameters()
{
  /// @brief Sets all simulation parameters from geo file

  //Parameters to be read in as strings
  std::string simStep="timeStep";
  std::string totNoFrames="totalNoFrames";

  std::string boundingBoxPos="gridOrigin";
  std::string boundingBoxSize="gridSize";
  std::string noCells="noGridCells";

  std::string lameMu="LameMu";
  std::string lameLambda="LameLambda";
  std::string compLimit="CompressionLimit";
  std::string stretchLimit="StretchLimit";
  std::string hardnessCoeff="HardnessCoefficient";

  std::string heatCapSolid="HeatCapacitySolid";
  std::string heatCapFluid="HeatCapacityFluid";
  std::string heatCondSolid="HeatConductivitySolid";
  std::string heatCondFluid="HeatConductivityFluid";
  std::string latentHeat="LatentHeat";
  std::string freezeTemp="FreezingTemperature";

  std::string ambientTemp="ambientTemperature";
  std::string heatSourceTemp="heatSourceTemperature";

  //Need to read in each value from file
  ReadGeo* file=new ReadGeo(m_readFileName);

  m_simTimeStep=file->getSimulationParameter_Float(simStep);
  m_totalNoFrames=file->getSimulationParameter_Float(totNoFrames);

  m_boundingBoxPosition=file->getSimulationParameter_Vec3(boundingBoxPos);
  m_boundingBoxSize=file->getSimulationParameter_Float(boundingBoxSize);
  m_noCells=file->getSimulationParameter_Float(noCells);

  m_lameMuConstant=file->getSimulationParameter_Float(lameMu);
  m_lameLambdaConstant=file->getSimulationParameter_Float(lameLambda);
  m_compressionLimit=file->getSimulationParameter_Float(compLimit);
  m_stretchLimit=file->getSimulationParameter_Float(stretchLimit);
  m_hardnessCoefficient=file->getSimulationParameter_Float(hardnessCoeff);

  m_heatCapacitySolid=file->getSimulationParameter_Float(heatCapSolid);
  m_heatCapacityFluid=file->getSimulationParameter_Float(heatCapFluid);
  m_heatConductivitySolid=file->getSimulationParameter_Float(heatCondSolid);
  m_heatConductivityFluid=file->getSimulationParameter_Float(heatCondFluid);
  m_latentHeat=file->getSimulationParameter_Float(latentHeat);
  m_freezingTemperature=file->getSimulationParameter_Float(freezeTemp);

  m_ambientTemperature=file->getSimulationParameter_Float(ambientTemp);
  m_heatSourceTemperature=file->getSimulationParameter_Float(heatSourceTemp);

  delete file;

}

//----------------------------------------------------------------------------------------------------------------------

void SimulationController::setupParticles()
{
  //Set up vectors to contain positions, mass, phase and temperature
  std::string mass="mass";
  std::string phase="phase";
  std::string temperature="temperature";

  std::vector<Eigen::Vector3f> positionList;
  std::vector<float> massList;
  std::vector<float> phaseList;
  std::vector<float> temperatureList;

  //Read in the data from file
  ReadGeo* file=new ReadGeo(m_readFileName);

  file->getPointPositions(m_noParticles, positionList);
  file->getPointParameter_Float(mass, massList);
  file->getPointParameter_Float(phase, phaseList);
  file->getPointParameter_Float(temperature, temperatureList);

  delete file;

  //Create emitter by passing in the data
  m_noParticles=positionList.size();
  m_emitter->createParticles(m_noParticles, positionList, massList, temperatureList, phaseList);


}

//----------------------------------------------------------------------------------------------------------------------

void SimulationController::update()
{
  /// @brief Steps the simulation. This controls the interlink between the particles and the grid

  //Determine if first step
  bool isFirstStep=false;
  if (m_noFrames==0 && m_elapsedTimeAfterFrame==0)
  {
    isFirstStep=true;
  }

  //Update elastic/plastic
  m_emitter->presetParticles(m_velocityContributionAlpha, m_temperatureContributionBeta);

  //Update grid which includes
  //Calculate interpolation weights
  //Transfer data from particles to grid
  //Calculate new velocity and temperature
  //Transfer data back to particles
  m_grid->update(m_simTimeStep, m_emitter, isFirstStep, m_velocityContributionAlpha, m_temperatureContributionBeta);


  //Update particles
  m_emitter->updateParticles(m_simTimeStep);

}

//----------------------------------------------------------------------------------------------------------------------

void SimulationController::render(ngl::Mat4 _modelMatrixCamera)
{
  /// @brief Renders particles through the emitter.Checks whether it is time to write to frame

  //Render particles
  m_emitter->renderParticles(_modelMatrixCamera, m_camera);

  //Check whether to create frame

}
