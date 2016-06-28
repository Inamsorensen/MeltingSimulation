#include "Emitter.h"

Emitter::Emitter(int _noParticles, float _particleMass)
{
  /// @brief Create particles. Need to set parameters separately

  m_noParticles=_noParticles;

  //Create particles
  for (int i=0; i<m_noParticles; i++)
  {
    Particle* particle=new Particle();
    particle->m_mass=_particleMass;
    m_particles.push_back(particle);
  }

}

//----------------------------------------------------------------------------------------------------------------------

Emitter::~Emitter()
{
  /// @brief Removes particle pointers

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
}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::setStrainConstants(float _lameMuConstant, float _lameLambdaConstant, float _compressionLim, float _stretchLim)
{
  /// @brief Set values to input values

  m_lameMuConstant=_lameMuConstant;
  m_lameLambdaConstant=_lameLambdaConstant;
  m_compressionLimit=_compressionLim;
  m_stretchLimit=_stretchLim;
}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::setTemperatureConstants(float _heatCapSolid, float _heatCapFluid, float _heatCondSolid, float _heatCondFluid, float _latentHeat, float _freezeTemp)
{
  /// @brief Set values to input values

  m_heatCapacitySolid=_heatCapSolid;
  m_heatCapacityFluid=_heatCapFluid;
  m_heatConductivitySolid=_heatCondSolid;
  m_heatConductivityFluid=_heatCondFluid;
  m_latentHeat=_latentHeat;
  m_freezingTemperature=_freezeTemp;
}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::setParticlePosition(std::vector<ngl::Vec3> *_particlePositions)
{
  if (_particlePositions->size()==((unsigned)m_noParticles))
  {
    for (int i=0; i<m_noParticles; i++)
    {
      m_particles[i]->m_position=_particlePositions->at(i);
    }
  }
  else
  {
    std::cout<<"Different number of particles in emitter and particle positions given.\n";
  }
}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::updateParticles()
{
  //Probably need to have more than one of these functions to say what's being updated. Also need transfer of data
}

//----------------------------------------------------------------------------------------------------------------------

void Emitter::renderParticles()
{
  /// @brief Render particles

  //For each particle, call render
  for (int i=0; i<m_noParticles; i++)
  {
    //Render particle
  }
}
