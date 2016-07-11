#include "Particle.h"

#include <math.h>

//----------------------------------------------------------------------------------------------------------------------

Particle::Particle(Eigen::Vector3f _position, float _mass, float _temperature, bool _isSolid, float _latentHeat, Emitter* _emitter)
{
  /// @brief Initiates particle values that need to be initialised

  m_position=_position;
  m_mass=_mass;
  m_temperature=_temperature;
  m_emitter=_emitter;


  //Initialise everything else to zero
  m_velocity.setZero();
  m_initialDensity=0.0;
  m_initialVolume=0.0;

  m_deformationElastic.setZero();
  m_deformationPlastic.setZero();
  m_deformationElastic_Deviatoric.setZero();
  m_deformationElastic_Deviatoric_Diff.setZero();
  m_R_deformationElastic_Deviatoric.setZero();
  m_S_deformationElastic_Deviatoric.setZero();


  //Set these to ones so don't divide by zero
  m_detDeformGrad=1.0;
  m_detDeformGradElastic=1.0;
  m_detDeformGradPlastic=1.0; //Not sure if should be one.

  //Set lame values to one for now since transfer 1/lameLambda. Updated at beginning of update step
  m_lameMu=1.0;
  m_lameLambda=1.0;


  //Set phase and fill/empty transition heat depending on whether solid or liquid
  if (_isSolid)
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

void Particle::getParticleData_CellFace(float &o_mass, Eigen::Vector3f &o_velocity, Phase &o_phase)
{
  /// @brief Used to collect data that is to be transferred to the grid cell faces

  o_mass=m_mass;
  o_velocity=m_velocity;
  o_phase=m_phase;

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::getParticleData_CellCentre(float &o_mass, float &o_detDeformGrad, float &o_detDeformGradElast, Phase &o_phase, float &o_temp, float &o_lameLambdaInverse)
{
  o_mass=m_mass;
  o_detDeformGrad=m_detDeformGrad;
  o_detDeformGradElast=m_detDeformGradElastic;
  o_phase=m_phase;
  o_temp=m_temperature;
  o_lameLambdaInverse=(1.0/m_lameLambda);

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
