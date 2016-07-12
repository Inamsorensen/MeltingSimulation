#include "Particle.h"

#include <math.h>
#include <algorithm>

#include "MathFunctions.h"
#include "Emitter.h"

//----------------------------------------------------------------------------------------------------------------------

Particle::Particle(Eigen::Vector3f _position, float _mass, float _temperature, bool _isSolid, float _latentHeat, Emitter* _emitter)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Sets position, mass, temperature, phase, latent heat and emitter pointer to the function input parameters

  For remaining class variables, set all to zero. Except determinants of deformation gradients, and
  lame coefficients since these will be used in division.

  All the values set to zero should be set using separate functions. Eg. setLameCoefficients

  Depending on phase the transition heat is either set to latent heat for liquid, and zero for solid.

  ------------------------------------------------------------------------------------------------------
  */

  m_position=_position;
  m_mass=_mass;
  m_temperature=_temperature;
  m_emitter=_emitter;


  //Initialise everything else to zero
  m_velocity.setZero();
  m_initialDensity=0.0;
  m_initialVolume=0.0;

  //Initialise matrices to identity matrix
  m_deformationElastic=m_deformationElastic.Identity();
  m_deformationPlastic=m_deformationPlastic.Identity();
  m_deformationElastic_Deviatoric=m_deformationElastic_Deviatoric.Identity();
  m_R_deformationElastic_Deviatoric=m_R_deformationElastic_Deviatoric.Identity();
  m_S_deformationElastic_Deviatoric=m_S_deformationElastic_Deviatoric.Identity();

  //Set these to ones so don't divide by zero
  m_detDeformGrad=1.0;
  m_detDeformGradElastic=1.0;
  m_detDeformGradPlastic=1.0; //Not sure if should be one.

  //Set lame values to one for now since transfer 1/lameLambda. Updated at beginning of update step
  m_lameMu=1.0;
  m_lameLambda=1.0;

  //Set dimension to 3 for now. Might be that this should be 1
  m_dimension=3;

  //Set phase and fill/empty transition heat depending on whether solid or liquid
  if (_isSolid)
  {
    m_phase=Phase::Solid;
    m_transitionHeat=0.0;
  }
  else
  {
    m_phase=Phase::Liquid;
    m_transitionHeat=_latentHeat;
  }

}

//----------------------------------------------------------------------------------------------------------------------

Particle::~Particle()
{

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::getParticleData_CellFace(float &o_mass, Eigen::Vector3f &o_velocity, Phase &o_phase)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Used to collect data to be transferred to the grid cell faces

  Improvements: Could have single functions for each to make more versatile, or make friend class
  ------------------------------------------------------------------------------------------------------
  */

  o_mass=m_mass;
  o_velocity=m_velocity;
  o_phase=m_phase;

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::getParticleData_CellCentre(float &o_mass, float &o_detDeformGrad, float &o_detDeformGradElast, Phase &o_phase, float &o_temp, float &o_lameLambdaInverse)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Used to collect data to be transferred to the grid cell centres

  Improvements: Could have single functions for each to make more versatile, or make friend class
  ------------------------------------------------------------------------------------------------------
  */

  o_mass=m_mass;
  o_detDeformGrad=m_detDeformGrad;
  o_detDeformGradElast=m_detDeformGradElastic;
  o_phase=m_phase;
  o_temp=m_temperature;
  o_lameLambdaInverse=(1.0/m_lameLambda);

}

//----------------------------------------------------------------------------------------------------------------------

Eigen::Matrix3f Particle::getDeformEDevDiff_Z(const Eigen::Matrix3f &_Z)
{
  /* Outline - Based on technical paper
  ------------------------------------------------------------------------------------------------------
  B:Z where B is the differential of J^{-1/d}F is given by

  JE^{-1/d}*[Z+(-1/d)(FE^{-T}:Z)FE]
  ------------------------------------------------------------------------------------------------------
  */

  //Set up result matrix
  Eigen::Matrix3f result;

  //Calculate -1/d
  float dimInverse=(-1.0/((float)m_dimension));

  //Calculate JE^{-1/d}
  float detDeformGrad_dimInv=pow(m_detDeformGradElastic, dimInverse);

  //Calculate FE^{-T}
  Eigen::Matrix3f deformGrad_Trans=m_deformationElastic.transpose();
  Eigen::Matrix3f deformGrad_TransInv=deformGrad_Trans.inverse();

  //Calculate FE^{-T}:Z
  Eigen::Matrix3f deformGrad_TransInv_Z=MathFunctions::matrixElementMultiplication(deformGrad_TransInv, _Z);

  //Calculate (FE^{-T}:Z)FE
  result=deformGrad_TransInv_Z*m_deformationElastic;

  //Calculate (-1/d)(FE^{-T}:Z)FE
  result*=dimInverse;

  //Add Z
  result+=_Z;

  //Multiply by JE^{-1/d}
  result*=detDeformGrad_dimInv;

  //Return result
  return result;

}

//----------------------------------------------------------------------------------------------------------------------

Eigen::Matrix3f Particle::getZ_DeformEDevDiff(const Eigen::Matrix3f &_Z)
{
  /* Outline - Based on technical paper
  ------------------------------------------------------------------------------------------------------
  Z:B where B is the differential of J^{-1/d}F is given by

  JE^{-1/d}*[Z+(-1/d)(FE:Z)FE^{-T}]
  ------------------------------------------------------------------------------------------------------
  */

  //Set up result matrix
  Eigen::Matrix3f result;

  //Calculate -1/d
  float dimInverse=(-1.0/((float)m_dimension));

  //Calculate JE^{-1/d}
  float detDeformGrad_dimInv=pow(m_detDeformGradElastic, dimInverse);

  //Calculate FE^{-T}
  Eigen::Matrix3f deformGrad_Trans=m_deformationElastic.transpose();
  Eigen::Matrix3f deformGrad_TransInv=deformGrad_Trans.inverse();

  //Calculate FE:Z
  Eigen::Matrix3f deformGrad_Z=MathFunctions::matrixElementMultiplication(m_deformationElastic, _Z);

  //Calculate (FE:Z)FE^{-1/d}
  result=deformGrad_Z*deformGrad_TransInv;

  //Calculate (-1/d)(FE^{-T}:Z)FE
  result*=dimInverse;

  //Add Z
  result+=_Z;

  //Multiply by JE^{-1/d}
  result*=detDeformGrad_dimInv;

  //Return result
  return result;

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::presetParticlesForTimeStep()
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Calculate plasticity contribution

  Apply correction for splitting
    FE=JP^{1/d}FE, FP=JP^{-1/d}FP

  Calculate determinants and store
  NB! Not sure whether determinants should be calculated before or after correction for splitting

  Calculate new lame coefficients

  Calculate deviatoric elastic matrices + polar decomposition
  ------------------------------------------------------------------------------------------------------
  */

  //Apply plasticity contribution
  applyPlasticity();

  //Apply correction for splitting
  m_detDeformGradPlastic=m_deformationPlastic.determinant();

  //For test set determinant to absolute value. Determinant should always be positive
  ///NB!!!!!!!! Not sure I can do this is reality
  if (m_detDeformGradPlastic<0)
  {
    std::cout<<"Determinant of F_{P} is negative. Setting it to positive.\n";
    m_detDeformGradPlastic=std::abs(m_detDeformGradPlastic);
  }

  float dimensionInv=1.0/((float)m_dimension);

  float plasticCorrectionForElastic=pow(m_detDeformGradPlastic, dimensionInv);
  float plasticCorrectionForPlastic=pow(m_detDeformGradPlastic, -dimensionInv);

  m_deformationElastic*=plasticCorrectionForElastic;
  m_deformationPlastic*=plasticCorrectionForPlastic;

  //Calculate new determinants - NB! Not sure if this should be done before splitting correction
  Eigen::Matrix3f deformationGradient=m_deformationElastic*m_deformationPlastic;
  m_detDeformGrad=deformationGradient.determinant();
  m_detDeformGradElastic=m_deformationElastic.determinant();
  m_detDeformGradPlastic=m_deformationPlastic.determinant();

  //For test set determinant to absolute value. Determinant should always be positive
  ///NB!!!!!!!! Not sure I can do this is reality
  if (m_detDeformGrad<0)
  {
    std::cout<<"Determinant of F is negative. Setting it to positive.\n";
    m_detDeformGrad=std::abs(m_detDeformGrad);
  }
  if (m_detDeformGradElastic<0)
  {
    std::cout<<"Determinant of F_{E} is negative. Setting it to positive.\n";
    m_detDeformGradElastic=std::abs(m_detDeformGradElastic);
  }
  if (m_detDeformGradPlastic<0)
  {
    std::cout<<"Determinant of F_{P} is negative. Setting it to positive.\n";
    m_detDeformGradPlastic=std::abs(m_detDeformGradPlastic);
  }


  //Calculate new lame coefficients
  float lameMuConstant=m_emitter->getLameMuConstant();
  float lameLambdaConstant=m_emitter->getLameLambdaConstant();
  float hardness=m_emitter->getHardnessCoefficient();
  float exponentialParam=hardness*(1.0-m_detDeformGradPlastic);
  float hardnessImpact=exp(exponentialParam);

  if (m_phase==Phase::Liquid)
  {
    m_lameMu=0.0;
  }
  else
  {
    m_lameMu=lameMuConstant*hardnessImpact;
  }
  m_lameLambda=lameLambdaConstant*hardnessImpact;


  //Calculate new elastic deviatoric components
  float elasticCorrectionForElastic=pow(m_detDeformGradElastic, -dimensionInv);
  m_deformationElastic_Deviatoric=elasticCorrectionForElastic*m_deformationElastic;

  MathFunctions::polarDecomposition(m_deformationElastic_Deviatoric, m_R_deformationElastic_Deviatoric, m_S_deformationElastic_Deviatoric);


}

//----------------------------------------------------------------------------------------------------------------------

void Particle::update(float _dt)
{

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::applyPlasticity()
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Calculate F^{n+1}
  Calculate singular value decomp of FE
  Clamp values to compression and stretch limits using std::min, std::max
  Recompute FE
  Calculate new FP
  Verify F^{n+1}==FE*FP
  ------------------------------------------------------------------------------------------------------
  */

  //Calculate F^{n+1}
  //deformationElastic has been updated here from end of previous step whilst deformationPlastic was updated
  //at the beginning of the previous step. Hence n+1 and n
  Eigen::Matrix3f deformationGradient=m_deformationElastic*m_deformationPlastic;

  //Singular value decomposition
  Eigen::Matrix3f U_matrix;
  Eigen::Matrix3f V_matrix;
  Eigen::Matrix3f U_matrix_trans;
  Eigen::Matrix3f V_matrix_trans;
  Eigen::Matrix3f singularValue_matrix;

  MathFunctions::singularValueDecomposition(m_deformationElastic, U_matrix, singularValue_matrix, V_matrix);

  U_matrix_trans=U_matrix.transpose();
  V_matrix_trans=V_matrix.transpose();


  //Get compression and stretch limits
  float compressionLimit=m_emitter->getCompressionLimit();
  float stretchLimit=m_emitter->getStretchLimit();


  //Clamp singular values
  for (int i=0; i<3; ++i)
  {
    float clampedValue=singularValue_matrix(i,i);

    //Check if greater than 1+stretch limit by using std::min. If is larger, then will return stretchlimit+1
    clampedValue=std::min<float>(clampedValue, (1.0+stretchLimit));

    //Check if smaller than 1-compression limit by using std::max. If smaller, then will return 1-compression limit
    clampedValue=std::max<float>(clampedValue, (1.0-compressionLimit));

    //Reinsert clamped value
    singularValue_matrix(i,i)=clampedValue;
  }


  //Recalculate FE
  Eigen::Matrix3f tempMatrix=singularValue_matrix*V_matrix_trans;
  m_deformationElastic=U_matrix*tempMatrix;


  //Calculate new FP
  Eigen::Matrix3f singularValue_matrix_inverse=singularValue_matrix.inverse();
  tempMatrix=V_matrix*singularValue_matrix_inverse;
  Eigen::Matrix3f tempMatrix2=tempMatrix*U_matrix_trans;
  m_deformationPlastic=tempMatrix2*deformationGradient;


  //Verify that new deformation gradient is the same as the one before applying plasticity
  Eigen::Matrix3f newDeformationGradient=m_deformationElastic*m_deformationPlastic;

  //Check if new deformation gradient is the same as the one calculated at the beginning
  //Even when correct there is a difference probably due to floating point error.
//  if (newDeformationGradient!=deformationGradient)
//  {
//    throw std::invalid_argument("New deformation gradient not the same as the old.");
//  }


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
