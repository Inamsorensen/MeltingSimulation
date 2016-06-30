#include "Particle.h"

//----------------------------------------------------------------------------------------------------------------------

Particle::Particle(Eigen::Vector3f _position, float _mass, float _temperature, bool _solid, float _latentHeat, Emitter *_emitter)
{
  /// @brief Initiates particle values that need to be initialised

  m_position=_position;
  m_mass=_mass;
  m_temperature=_temperature;
  m_emitter=_emitter;

  //Set phase and fill/empty transition heat depending on whether solid or liquid
  if (_solid==true)
  {
    m_phase=Phase::Solid;
    m_transitionHeat=_latentHeat;
  }
  else
  {
    m_phase=Phase::Liquid;
    m_transitionHeat=0.0;
  }

}

//----------------------------------------------------------------------------------------------------------------------

Particle::~Particle()
{

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::update(float _dt)
{

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::updateVelocity()
{

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::updatePosition()
{

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::updateTemperature()
{

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::updateDeformationGradient()
{

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::render()
{

}

//----------------------------------------------------------------------------------------------------------------------
