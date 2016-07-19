#include <ngl/Transformation.h>
#include <ngl/ShaderLib.h>
#include <ngl/VAOPrimitives.h>

#include "Emitter.h"

Emitter::Emitter()
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Initialise all values to zero

  Set all of these variables separately using
    createParticles
    setStrainConstants
    setRenderParameters

  ------------------------------------------------------------------------------------------------------
  */

  m_noParticles=0;

  m_lameMuConstant=0.0;
  m_lameLambdaConstant=0.0;
  m_hardnessCoefficient=0.0;
  m_compressionLimit=0.0;
  m_stretchLimit=0.0;
  m_heatCapacitySolid=0.0;
  m_heatCapacityFluid=0.0;
  m_heatConductivitySolid=0.0;
  m_heatConductivityFluid=0.0;
  m_latentHeat=0.0;
  m_freezingTemperature=0.0;
  m_freezingTemperatureBuffer=2.0;

  m_xMin=0.0;
  m_xMax=0.0;
  m_yMin=0.0;
  m_yMax=0.0;
  m_zMin=0.0;
  m_zMax=0.0;

  m_particleShaderName="";
  m_particleRadius=0.0;

}

//----------------------------------------------------------------------------------------------------------------------

Emitter::~Emitter()
{
  /* Outline
  ------------------------------------------------------------------------------------------------------

  Deletes particle pointers

  ------------------------------------------------------------------------------------------------------
  */

  //Check number of particles in vector
  int noParticlesCurrent=m_particles.size();

  //If more than zero particle pointers in vector, delete these pointers
  if (noParticlesCurrent!=0)
  {
    std::cout<<"Deleting particles\n";

    for (int i=0; i<noParticlesCurrent; i++)
    {
      delete m_particles[i];
    }
  }

  //Clear vector
  m_particles.clear();

  std::cout<<"Deleting emitter\n";
}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::createParticles(int _noParticles, const std::vector<Eigen::Vector3f> &_particlePositions, const std::vector<float> &_particleMass, const std::vector<float> &_particleTemperature, const std::vector<float> &_particlePhase)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Generates particles based on input parameters.
  ------------------------------------------------------------------------------------------------------
  */

  m_noParticles=_noParticles;

  //Create particles
  for (int i=0; i<m_noParticles; i++)
  {
    Eigen::Vector3f position=_particlePositions.at(i);
    float mass=_particleMass.at(i);
    float temperature=_particleTemperature.at(i)+273.0;  //Add 273 as temperature in Kelvin whereas read in is in Celsius
    bool solid=_particlePhase.at(i);

    Particle* particle=new Particle(position, mass, temperature, solid, m_latentHeat, this);
    m_particles.push_back(particle);
  }
}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::setStrainConstants(float _lameMuConstant, float _lameLambdaConstant, float _compressionLim, float _stretchLim, float _hardnessCoefficient)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Sets strain constants from input
  ------------------------------------------------------------------------------------------------------
  */

  m_lameMuConstant=_lameMuConstant;
  m_lameLambdaConstant=_lameLambdaConstant;
  m_compressionLimit=_compressionLim;
  m_stretchLimit=_stretchLim;
  m_hardnessCoefficient=_hardnessCoefficient;

}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::setTemperatureConstants(float _heatCapSolid, float _heatCapFluid, float _heatCondSolid, float _heatCondFluid, float _latentHeat, float _freezeTemp)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Set temperature constants from input
  ------------------------------------------------------------------------------------------------------
  */

  m_heatCapacitySolid=_heatCapSolid;
  m_heatCapacityFluid=_heatCapFluid;
  m_heatConductivitySolid=_heatCondSolid;
  m_heatConductivityFluid=_heatCondFluid;
  m_latentHeat=_latentHeat;
  m_freezingTemperature=_freezeTemp+273.0; //Add 273 to go from Celsius to Kelvin
}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::setCollisionObject(float _xMin, float _xMax, float _yMin, float _yMax, float _zMin, float _zMax)
{
  m_xMin=_xMin;
  m_xMax=_xMax;
  m_yMin=_yMin;
  m_yMax=_yMax;
  m_zMin=_zMin;
  m_zMax=_zMax;
}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::setRenderParameters(std::string _shaderName, float _particleRadius)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Set render parameters from input
  ------------------------------------------------------------------------------------------------------
  */

  m_particleShaderName=_shaderName;
  m_particleRadius=_particleRadius;
}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::presetParticles(float _velocityContribAlpha, float _tempContribBeta)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Calls particle preset for all particles
  ------------------------------------------------------------------------------------------------------
  */

#pragma omp parallel for
  for (int i=0; i<m_noParticles; ++i)
  {
    m_particles[i]->presetParticlesForTimeStep(_velocityContribAlpha, _tempContribBeta);
  }
}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::updateParticles(float _dt)
{
#pragma omp parallel for
  for (int i=0; i<m_noParticles; i++)
  {
    m_particles[i]->update(_dt, m_xMin, m_xMax, m_yMin, m_yMax, m_zMin, m_zMax);
  }
}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::renderParticles(ngl::Mat4 _modelMatrixCamera, ngl::Camera* _camera)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Render particles contained by emitter
  ------------------------------------------------------------------------------------------------------
  */

  //Get shader and VAO
  ngl::ShaderLib* shaderLib=ngl::ShaderLib::instance();
  ngl::VAOPrimitives* vaoPrimitives=ngl::VAOPrimitives::instance();

  //Set shader to use
  shaderLib->use(m_particleShaderName);


  for (int i=0; i<m_noParticles; i++)
  {
    //Get particle position and make into Vec4
    ngl::Vec4 particlePosition;
    Eigen::Vector3f particlePositionVec3=m_particles[i]->getPosition();
    particlePosition.m_x=particlePositionVec3(0);
    particlePosition.m_y=particlePositionVec3(1);
    particlePosition.m_z=particlePositionVec3(2);
    particlePosition.m_w=1.0;

    //Calculate new position due to scene transformations
    ngl::Vec4 transformedPosition=particlePosition*_modelMatrixCamera;

    //Set transformation matrix
    ngl::Transformation transformationMatrix;
    transformationMatrix.setPosition(transformedPosition.m_x, transformedPosition.m_y, transformedPosition.m_z);
    transformationMatrix.setScale(m_particleRadius, m_particleRadius, m_particleRadius);

    //Calculate MVP matrices
    ngl::Mat4 M;
    ngl::Mat4 MV;
    ngl::Mat4 MVP;
    ngl::Mat3 normalMatrix;

    M=transformationMatrix.getMatrix();
    MV=M*_camera->getViewMatrix();
    MVP=MV*_camera->getProjectionMatrix();
    normalMatrix=MV;
    normalMatrix.inverse();

    //Set MVP to shader
    shaderLib->setShaderParamFromMat4("M", M);
    shaderLib->setShaderParamFromMat4("MV", MV);
    shaderLib->setShaderParamFromMat4("MVP", MVP);
    shaderLib->setShaderParamFromMat3("normalMatrix", normalMatrix);

    //Draw particle
    vaoPrimitives->draw("sphere");
  }
}
