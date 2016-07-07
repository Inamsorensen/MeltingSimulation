#include  <iostream>

#include "SimulationController.h"

//----------------------------------------------------------------------------------------------------------------------

SimulationController* SimulationController::m_instance=nullptr;

//----------------------------------------------------------------------------------------------------------------------

SimulationController::SimulationController()
{
  /// @brief Set up simulation parameters and create grid and emitter pointers

  //Setup if don't read from file
  //Timer setup
  m_simTimeStep=0.01;
  m_elapsedTimeAfterFrame=0.0;
  m_noFrames=0;

  //Grid setup
  Eigen::Vector3f gridPos;
  gridPos(0)=-0.5;
  gridPos(1)=-0.5;
  gridPos(2)=-0.5;
  m_gridPosition=gridPos;
  m_gridSize=1.0;
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

  //Read in simulation parameters
  std::string simulationParametersFile="../HoudiniFiles/particles.geo";
  readSimulationParameters(simulationParametersFile);

  //Create emitter and particles
  m_emitter=new Emitter();
  m_emitter->setStrainConstants(m_lameMuConstant, m_lameLambdaConstant, m_compressionLimit, m_stretchLimit);
  m_emitter->setTemperatureConstants(m_heatCapacitySolid, m_heatCapacityFluid, m_heatConductivitySolid, m_heatConductivityFluid, m_latentHeat, m_freezingTemperature);
  setupParticles();

  //Create grid
  //Need to stagger grid as Houdini setup has origin in lower back corner, but MAC staggered
  //has origin in the middle of the cell in the lower back corner
  float halfCellSize=(1.0/2.0)*(m_gridSize/((float)m_noCells));
  Eigen::Vector3f staggeredGridPosition=m_gridPosition;
  staggeredGridPosition(0)+=halfCellSize;
  staggeredGridPosition(1)+=halfCellSize;
  staggeredGridPosition(2)+=halfCellSize;
  m_grid=Grid::createGrid(staggeredGridPosition, m_gridSize, m_noCells);
  m_grid->findParticleInCell(m_emitter);
  m_grid->TEST_findParticleInCell(m_emitter);

  //std::cout<<"test\n";

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
  m_particleRadius=0.1;

  m_emitter->setRenderParameters(m_shaderName, m_particleRadius);

}

//----------------------------------------------------------------------------------------------------------------------

void SimulationController::readSimulationParameters(std::string _fileName)
{
  /// @brief Sets all simulation parameters from geo file

  //Parameters to be read in as strings
  std::string simStep="timeStep";
  std::string totNoFrames="totalNoFrames";

  std::string gridPos="gridOrigin";
  std::string gridSize="gridSize";
  std::string noCells="noGridCells";

  std::string lameMu="LameMu";
  std::string lameLambda="LameLambda";
  std::string compLimit="CompressionLimit";
  std::string stretchLimit="StretchLimit";

  std::string heatCapSolid="HeatCapacitySolid";
  std::string heatCapFluid="HeatCapacityFluid";
  std::string heatCondSolid="HeatConductivitySolid";
  std::string heatCondFluid="HeatConductivityFluid";
  std::string latentHeat="LatentHeat";
  std::string freezeTemp="FreezingTemperature";

  //Need to read in each value from file
  ReadGeo* file=new ReadGeo(_fileName);

  m_simTimeStep=file->getSimulationParameter_Float(simStep);
  m_totalNoFrames=file->getSimulationParameter_Float(totNoFrames);

  m_gridPosition=file->getSimulationParameter_Vec3(gridPos);
  m_gridSize=file->getSimulationParameter_Float(gridSize);
  m_noCells=file->getSimulationParameter_Float(noCells);

  m_lameMuConstant=file->getSimulationParameter_Float(lameMu);
  m_lameLambdaConstant=file->getSimulationParameter_Float(lameLambda);
  m_compressionLimit=file->getSimulationParameter_Float(compLimit);
  m_stretchLimit=file->getSimulationParameter_Float(stretchLimit);

  m_heatCapacitySolid=file->getSimulationParameter_Float(heatCapSolid);
  m_heatCapacityFluid=file->getSimulationParameter_Float(heatCapFluid);
  m_heatConductivitySolid=file->getSimulationParameter_Float(heatCondSolid);
  m_heatConductivityFluid=file->getSimulationParameter_Float(heatCondFluid);
  m_latentHeat=file->getSimulationParameter_Float(latentHeat);
  m_freezingTemperature=file->getSimulationParameter_Float(freezeTemp);

  delete file;

}

//----------------------------------------------------------------------------------------------------------------------

void SimulationController::setupParticles()
{
  std::string particleFileName="../HoudiniFiles/particles.geo";

  //Set up vectors to contain positions, mass, phase and temperature
  std::string mass="mass";
  std::string phase="phase";
  std::string temperature="temperature";

  std::vector<Eigen::Vector3f> positionList;
  std::vector<float> massList;
  std::vector<float> phaseList;
  std::vector<float> temperatureList;

  //Read in the data from file
  ReadGeo* file=new ReadGeo(particleFileName);

  file->getPointPositions(&positionList);
  file->getPointParameter_Float(mass, &massList);
  file->getPointParameter_Float(phase, &phaseList);
  file->getPointParameter_Float(temperature, &temperatureList);

  delete file;

  //Create emitter by passing in the data
  int noParticles=positionList.size();
  m_emitter->createParticles(noParticles, &positionList, &massList, &temperatureList, &phaseList);


}

//----------------------------------------------------------------------------------------------------------------------

void SimulationController::update()
{
  /// @brief Steps the simulation. This controls the interlink between the particles and the grid

  //Update elastic/plastic

  //Calculate interpolation weights

  //Transfer data from particles to grid

  //Update grid to calculate new velocity and temperature

  //Transfer data back to particles

  //Update particles

}

//----------------------------------------------------------------------------------------------------------------------

void SimulationController::render(ngl::Mat4 _modelMatrixCamera)
{
  /// @brief Renders particles through the emitter.Checks whether it is time to write to frame

  //Render particles
  m_emitter->renderParticles(_modelMatrixCamera, m_camera);

  //Check whether to create frame

}
