#include "Particle.h"

#include <math.h>
#include <algorithm>
#include <cmath>

#include "MathFunctions.h"
#include "Emitter.h"

//----------------------------------------------------------------------------------------------------------------------

Particle::Particle(unsigned int _id, Eigen::Vector3f _position, float _mass, float _temperature, bool _isSolid, float _latentHeat, Emitter* _emitter)
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

  m_id=_id;

  m_position=_position;
  m_mass=_mass;
  m_temperature=_temperature;
  m_previousTemperature=_temperature;
  m_emitter=_emitter;


  //Initialise everything else to zero
  m_velocity.setZero();
  m_initialDensity=0.0;
  m_initialVolume=0.0;
  m_previousVelocity.setZero();
  m_velocityGradient=m_velocityGradient.Identity();

  //Initialise matrices to identity matrix
  m_deformationElastic=m_deformationElastic.Identity();
  m_deformationPlastic=m_deformationPlastic.Identity();
  m_deformationElastic_Deviatoric=m_deformationElastic_Deviatoric.Identity();
  m_R_deformationElastic_Deviatoric=m_R_deformationElastic_Deviatoric.Identity();
  m_S_deformationElastic_Deviatoric=m_S_deformationElastic_Deviatoric.Identity();

//  //NB!!!!!!!!!!!!!!!!!!!!!
//  //FOR TESTING APPLY PLASTICITY CALC ONLY
//  m_deformationElastic(0,0)=1.0;
//  m_deformationElastic(0,1)=2.0;
//  m_deformationElastic(0,2)=0.0;
//  m_deformationElastic(1,0)=2.0;
//  m_deformationElastic(1,1)=3.0;
//  m_deformationElastic(1,2)=3.0;
//  m_deformationElastic(2,0)=2.0;
//  m_deformationElastic(2,1)=2.0;
//  m_deformationElastic(2,2)=1.0;

//  m_deformationPlastic(0,0)=1.0;
//  m_deformationPlastic(0,1)=2.0;
//  m_deformationPlastic(0,2)=0.0;
//  m_deformationPlastic(1,0)=1.0;
//  m_deformationPlastic(1,1)=2.0;
//  m_deformationPlastic(1,2)=3.0;
//  m_deformationPlastic(2,0)=1.0;
//  m_deformationPlastic(2,1)=2.0;
//  m_deformationPlastic(2,2)=3.0;

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
  o_velocity=m_previousVelocity;
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
  o_temp=m_previousTemperature;
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

void Particle::presetParticlesForTimeStep(float _velocityContribAlpha, float _tempContribBeta)
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

  Calculate differentiated elasto-plastic energy
  ------------------------------------------------------------------------------------------------------
  */

  //Store old velocity and temperature.
  m_previousVelocity=m_velocity;
  m_previousTemperature=m_temperature;
  //Reset current velocity and temperature to alpha*vn and beta*Tn from FLIP update
  m_velocity=_velocityContribAlpha*m_previousVelocity;
  m_temperature=_tempContribBeta*m_previousTemperature;
  m_velocityGradient=m_velocityGradient.setZero();

  //Set new position to zero
  m_newPosition.setZero();


  //For now, update elastic deformation gradient if liquid, here
  if (m_phase==Phase::Liquid)
  {
    m_detDeformGradElastic=m_deformationElastic.determinant();
    float fluidCorrection=pow(m_detDeformGradElastic,(1.0/m_dimension));
    m_deformationElastic=m_deformationElastic.Identity();
    m_deformationElastic*=fluidCorrection;
  }


  //Apply plasticity contribution
  applyPlasticity();


  //For test set determinant to absolute value. Determinant should always be positive
//  ///NB!!!!!!!! Not sure I can do this is reality
//  if (m_detDeformGradPlastic<0)
//  {
//    std::cout<<"Determinant of F_{P} is negative. Setting it to positive.\n";
//    m_detDeformGradPlastic=std::abs(m_detDeformGradPlastic);
//  }


  //Calculate new lame coefficients
  float lameMuConstant=m_emitter->getLameMuConstant();
  float lameLambdaConstant=m_emitter->getLameLambdaConstant();
  float hardness=m_emitter->getHardnessCoefficient();
  float exponentialParam=hardness*(1.0-m_detDeformGradPlastic);
  float hardnessImpact=exp(exponentialParam);

  //Clamp hardness impact
  ///Not sure what to clamp to
  hardnessImpact=std::min<float>(hardnessImpact, 10.0);
  hardnessImpact=std::max<float>(hardnessImpact, 0.1);


  if (m_phase==Phase::Liquid)
  {
    m_lameMu=0.0; /// Unsure about this or whether I should just set FE=JE^(1/d)I? -> Should do both I believe
  }
  else
  {
    m_lameMu=lameMuConstant*hardnessImpact;
  }
  m_lameLambda=lameLambdaConstant*hardnessImpact;


  //Apply correction for splitting

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


//  //For test set determinant to absolute value. Determinant should always be positive
//  ///NB!!!!!!!! Not sure I can do this is reality
//  if (m_detDeformGrad<0)
//  {
//    std::cout<<"Determinant of F is negative. Setting it to positive.\n";
//    m_detDeformGrad=std::abs(m_detDeformGrad);
//  }
//  if (m_detDeformGradElastic<0)
//  {
//    std::cout<<"Determinant of F_{E} is negative. Setting it to positive.\n";
//    m_detDeformGradElastic=std::abs(m_detDeformGradElastic);
//  }
//  if (m_detDeformGradPlastic<0)
//  {
//    std::cout<<"Determinant of F_{P} is negative. Setting it to positive.\n";
//    m_detDeformGradPlastic=std::abs(m_detDeformGradPlastic);
//  }



//  //Calculate new lame coefficients
//  float lameMuConstant=m_emitter->getLameMuConstant();
//  float lameLambdaConstant=m_emitter->getLameLambdaConstant();
//  float hardness=m_emitter->getHardnessCoefficient();
//  float exponentialParam=hardness*(1.0-m_detDeformGradPlastic);
//  float hardnessImpact=exp(exponentialParam);

//  if (m_phase==Phase::Liquid)
//  {
//    m_lameMu=0.0; /// Unsure about this or whether I should just set FE=JE^(1/d)I? -> Should do both I believe
//  }
//  else
//  {
//    m_lameMu=lameMuConstant*hardnessImpact;
//  }
//  m_lameLambda=lameLambdaConstant*hardnessImpact;


  //Calculate new elastic deviatoric components
  float elasticCorrectionForElastic=pow(m_detDeformGradElastic, -dimensionInv);
  m_deformationElastic_Deviatoric=elasticCorrectionForElastic*m_deformationElastic;

  MathFunctions::polarDecomposition(m_deformationElastic_Deviatoric, m_R_deformationElastic_Deviatoric, m_S_deformationElastic_Deviatoric);


  //Calculate differential of elasto-plastic potential energy
  calcPotentialEnergyDiff();

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::update(float _dt, float _xMin, float _xMax, float _yMin, float _yMax, float _zMin, float _zMax)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Velocity and temperature is automatically updated from the grid

  Need to update deformation gradient

  Need to apply phase transitions. Ie. if temperature is above freezing temp but phase is solid and vice versa, then
  need transition stage

  Resolve collision between particles and surrounding objects

  Update particle positions
  ------------------------------------------------------------------------------------------------------
  */

  updateDeformationGradient(_dt);

  applyPhaseTransition();

  collisionResolve(_dt, _xMin, _xMax, _yMin, _yMax, _zMin, _zMax);

  updatePosition(_dt);
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


  //Calculate new determinants
  m_detDeformGrad=deformationGradient.determinant();
  m_detDeformGradElastic=m_deformationElastic.determinant();
  m_detDeformGradPlastic=m_deformationPlastic.determinant();

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::calcPotentialEnergyDiff()
{
  /* Outline
  ------------------------------------------------------------------------------------------------------

  dY_hatdFE=dYdFE_{p}:B_{kmij}

  dYdFE_{p}=2*lame_mu_{p}*JE^a*FE_{p} - 2*lame_mu_{p}*RE_{p} where a=-1/d where d=dimensions=3? for dimensions and RE is the polar decomposition of JE^a*FE

  B_{kmij}=JE^a*I + a*JE^a*FE_{p}^-T*FE_{p}

  ------------------------------------------------------------------------------------------------------
  */

  //Calculate dYdFE=2*mu*(FE-RE)
  Eigen::Matrix3f dYdFE=(2.0*m_lameMu)*(m_deformationElastic_Deviatoric-m_R_deformationElastic_Deviatoric);

  //Pass in dYdFE to particle's getZ_DeformEDevDiff. Return is dY^{hat}dFE
  m_potentialEnergyDiff=getZ_DeformEDevDiff(dYdFE);

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::updateDeformationGradient(float _dt)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  FE^{n+1}=R(dt*velGrad)FE^{n}

  velGrad (velocity gradient) is updated automatically by grid

  Need to determine R(dt*velGrad). Done in for loop to ensure that det(FE^{n+1})>0
  R(M)=I+M if det(I+M)>0 otherwise R(M)=R(M/2)^2
  ------------------------------------------------------------------------------------------------------
  */

  //Set first M and determinant
  Eigen::Matrix3f M_matrix=_dt*m_velocityGradient;

  //Identity matrix
  Eigen::Matrix3f I_matrix;
  I_matrix=I_matrix.Identity();

  //Determinant
  Eigen::Matrix3f determinantMatrix=I_matrix+M_matrix;
  float determinant=determinantMatrix.determinant();

  //Power count
  int powerCount=1;

  //Denominator
  float denom=(float(powerCount));

  //Stop loop
  bool stopLoop=false;

  if (determinant>0)
  {
    stopLoop=true;
  }

  while (!stopLoop)
  {
    //Increase power count
    powerCount+=1;

    //Calculate denominator for exponential series, ie. powerCount!
    denom*=(float(powerCount));

    //Calculate next component in exponential series
    Eigen::Matrix3f nextComponent=I_matrix;
    nextComponent*=(1.0/denom);

    for (int i=0; i<powerCount; i++)
    {
      nextComponent*=M_matrix;
    }

    //Calculate new determinant
    determinantMatrix+=nextComponent;
    determinant=determinantMatrix.determinant();

    //Check if determinant is positive
    if (determinant>0)
    {
      stopLoop=true;
    }
  }

  //Find new deformation gradient
  m_deformationElastic=determinantMatrix*m_deformationElastic;

}

//----------------------------------------------------------------------------------------------------------------------

void Particle::applyPhaseTransition()
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Get freezing temperature and latent heat and heat capacity

  Check if mismatch between new temperature and transition heat
    If T^{n+1}>T_freeze && L_p<L_trans
    if T^{n+1}<T_freeze && L_p>0

    If there is, increase transition heat by c*m*deltaT. And set temperature to freeze temp

    Check that L_p has not gone under 0 or over L_trans, if so clamp it
  ------------------------------------------------------------------------------------------------------
  */

  //Get transition temperature
  float transitionTemp=m_emitter->getTransitionTemperature();

  //Get latent heat
  float latentHeat=m_emitter->getLatentHeat();

  //Get heat capacity
  float heatCapacity;
  if (m_phase==Phase::Liquid)
  {
    heatCapacity=m_emitter->getHeatCapacityFluid();
  }
  else
  {
    heatCapacity=m_emitter->getHeatCapacitySolid();
  }

  //Check if previous temperature is at transition temp
  if (m_previousTemperature==transitionTemp)
  {
    //Increment transition heat
    float deltaT=m_temperature-m_previousTemperature;
    float heatIncrement=(heatCapacity*m_mass)*deltaT;
    m_transitionHeat+=heatIncrement;

    //Clamp transition heat
    m_transitionHeat=std::max<float>(m_transitionHeat, 0.0);
    m_transitionHeat=std::min<float>(m_transitionHeat, latentHeat);

    //Check if phase has changed
    if (m_transitionHeat==latentHeat)
    {
      m_phase=Phase::Liquid;
    }
    else if (m_transitionHeat==0.0)
    {
      m_phase=Phase::Solid;
    }
    else
    {
      //If transition heat has not been filled or emptied entirely, then no change in temperature
      m_temperature=transitionTemp;
    }
  }

  //Check if phase transition is reached
  else if ((m_previousTemperature<transitionTemp && m_temperature>=transitionTemp)
            || (m_previousTemperature>transitionTemp && m_temperature<=transitionTemp))
  {
    //Set temperature to transition temp, ignoring its contribution to transition heat
    m_temperature=transitionTemp;
  }
}

//----------------------------------------------------------------------------------------------------------------------

void Particle::collisionResolve(float _dt, float _xMin, float _xMax, float _yMin, float _yMax, float _zMin, float _zMax)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Check position of particle against all faces of bounding box
  ------------------------------------------------------------------------------------------------------
  */

  //Calculate possible new position
  Eigen::Vector3f possibleNewPosition=m_position+(_dt*m_velocity);
//  Eigen::Vector3f possibleNewPosition=m_newPosition;

  if (possibleNewPosition(0)<=_xMin || possibleNewPosition(0)>=_xMax || possibleNewPosition(1)<=_yMin || possibleNewPosition(1)>=_yMax || possibleNewPosition(2)<=_zMin || possibleNewPosition(2)>=_zMax)
  {
    m_velocity.setZero();
  }


}

//----------------------------------------------------------------------------------------------------------------------

void Particle::updatePosition(float _dt)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Update position using
  p^{n+1}=p^{n}+v*dt
  ------------------------------------------------------------------------------------------------------
  */

  m_position+=(_dt*m_velocity);
//  m_position=m_newPosition;
}

//----------------------------------------------------------------------------------------------------------------------
