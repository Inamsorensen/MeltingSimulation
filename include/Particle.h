#ifndef PARTICLE
#define PARTICLE

#include <eigen3/Eigen/Core>

#include <ngl/Vec3.h>
#include <ngl/Mat3.h>

//------------------------------------------------------------------------------------------------------------------------------------------------------
/// @file Particle.h
/// @brief Particle structure containing data specific to one particle
/// @author Ina M. Sorensen
/// @version 1.0
/// @date 25.06.16
///
/// @todo
//------------------------------------------------------------------------------------------------------------------------------------------------------

enum Phase
{
  Solid,
  Liquid
};

class Emitter;

class Particle
{
public:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Particle constructor
  //----------------------------------------------------------------------------------------------------------------------
  Particle(Eigen::Vector3f _position, float _mass, float _temperature, bool _isSolid, float _latentHeat, Emitter* _emitter);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Particle destructor
  //----------------------------------------------------------------------------------------------------------------------
  ~Particle();

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set Lame coefficients
  //----------------------------------------------------------------------------------------------------------------------
  void setLameCoefficients(float _lameMuConstant, float _lameLambdaConstant, float _hardnessCoefficient);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get particle position
  //----------------------------------------------------------------------------------------------------------------------
  inline Eigen::Vector3f getPosition(){return m_position;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get particle data for grid cell face
  //----------------------------------------------------------------------------------------------------------------------
  void getParticleData_CellFace(float &o_mass, Eigen::Vector3f &o_velocity, Phase &o_phase);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get particle data for grid cell centre
  //----------------------------------------------------------------------------------------------------------------------
  void getParticleData_CellCentre(float &o_mass, float &o_detDeformGrad, float &o_detDeformGradElast, Phase &o_phase, float &o_temp, float &o_lameLambdaInverse);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get particle data for grid cell centre
  //----------------------------------------------------------------------------------------------------------------------
  inline void addParticleDensity(float _densityIncrease){m_initialDensity+=_densityIncrease;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate initial volume
  //----------------------------------------------------------------------------------------------------------------------
  inline void calcInitialVolume(){m_initialVolume=m_mass/m_initialDensity;}

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Update particle. Calls to update velocity, position, temperature and deformation gradient
  /// @param [in] _dt: Time step
  //----------------------------------------------------------------------------------------------------------------------
  void update(float _dt);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Render particle
  //----------------------------------------------------------------------------------------------------------------------
  void render();

private:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Particle position
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Vector3f m_position;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Particle velocity
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Vector3f m_velocity;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Particle mass
  //----------------------------------------------------------------------------------------------------------------------
  float m_mass;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Initial density of particle. Used to calculate volume
  //----------------------------------------------------------------------------------------------------------------------
  float m_initialDensity;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Initial volume of particle.
  //----------------------------------------------------------------------------------------------------------------------
  float m_initialVolume;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Elastic deformation gradient, F_E
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Matrix3f m_deformationElastic;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Plastic deformation gradient, F_P;
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Matrix3f m_deformationPlastic;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Actual lame constant mu, taking into account hardening
  //----------------------------------------------------------------------------------------------------------------------
  float m_lameMu;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Actual lame constant lambda, taking into account hardening
  //----------------------------------------------------------------------------------------------------------------------
  float m_lameLambda;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Determinant of deformation gradient F
  //----------------------------------------------------------------------------------------------------------------------
  float m_detDeformGrad;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Determinant of elastic deformation gradient FE
  //----------------------------------------------------------------------------------------------------------------------
  float m_detDeformGradElastic;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Determinant of plastic deformation gradient FP
  //----------------------------------------------------------------------------------------------------------------------
  float m_detDeformGradPlastic;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief JE^(-1/d)*FE
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Matrix3f m_deformationElastic_Deviatoric;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief JE^(-1/d)*FE differentiated
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Matrix3f m_deformationElastic_Deviatoric_Diff;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief R of Polar decomposition of defElastic_Deviatoric
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Matrix3f m_R_deformationElastic_Deviatoric;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief S of Polar decomposition of defElastic_Deviatoric
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Matrix3f m_S_deformationElastic_Deviatoric;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Particle temperature
  //----------------------------------------------------------------------------------------------------------------------
  float m_temperature;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Transition heat of particle. Empty (0.0) if solid, full (equal to latent heat) if fluid. Transitioning if
  /// inbetween.
  //----------------------------------------------------------------------------------------------------------------------
  float m_transitionHeat;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Says if particle is solid or liquid.
  //----------------------------------------------------------------------------------------------------------------------
  Phase m_phase;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Emitter that particle belongs to
  //----------------------------------------------------------------------------------------------------------------------
  const Emitter* m_emitter;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Updates velocity
  //----------------------------------------------------------------------------------------------------------------------
  void updateVelocity();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Updates position
  //----------------------------------------------------------------------------------------------------------------------
  void updatePosition();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Updates temperature. Checks whether phase change is happening
  //----------------------------------------------------------------------------------------------------------------------
  void updateTemperature();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Updates deformation gradient and verifies elastic/plastic contribution
  //----------------------------------------------------------------------------------------------------------------------
  void updateDeformationGradient();

};


#endif // PARTICLE

