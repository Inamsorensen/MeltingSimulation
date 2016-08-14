#include "Grid.h"

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcDeviatoricForce_New(Particle* _particle, Eigen::Vector3f _eVector, Eigen::Vector3f _weightDiff)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  f_{i}=-sum_p(V_{p}*e_{a(i)}^T*dYdFE_{p}:B_{kmij}*FE_{p}^T*cubicBSpline_Diff_{ip}

  dYdFE_{p}=2*lame_mu_{p}*JE^a*FE_{p} - 2*lame_mu_{p}*RE_{p} where a=-1/d where d=3? for dimensions and RE is the polar decomposition of JE^a*FE

  B_{kmij}=JE^a*I + a*JE^a*FE_{p}^-T*FE_{p} - Stored in particle


  Calc dYdFE_{p}

  Special multiply dYdFE_{p} and deformationElastic_Deviatoric_Diff

  Multiply with FE^T

  Multiply with V_{p}

  Multiply with -1

  Return value

  ----------------------------------------------------------------------------------------------------------------
  */

  //Set final result
  float finalResult;

  //Get deformation gradient elastic, FE, and transpose it
  Eigen::Matrix3f deformationElastic=_particle->getDeformationElastic();
  Eigen::Matrix3f deformationElastic_trans=deformationElastic.transpose();

  //Get particle volume
  float particleVolume=_particle->getVolume();

  //Multiply FE^{T} with differentiated weight
  Eigen::Vector3f result_vec1=deformationElastic_trans*_weightDiff;

  //Get dY^{hat}dFE from particle
  Eigen::Matrix3f dY_hat_dFE=_particle->getPotentialEnergyDiff();

  //Multiply result with dY^{hat}dFE
  Eigen::Vector3f result_vec2=dY_hat_dFE*result_vec1;

  //Multiply result with e_{a(i)}
  finalResult=_eVector.dot(result_vec2);

  //Multiply result with particle volume
  finalResult*=particleVolume;

  //Make negative
  finalResult*=(-1.0);

  //Return result
  return finalResult;
}

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcBComponent_DeviatoricVelocity_New(Particle* _particle, Eigen::Vector3f _eVector, float _weight, float _deviatoricForce, float _massFace)
{
  /* Outline
  --------------------------------------------------------------------------------------------------------------
  b_{i}=w_{ip}m_{p}v_{p} + dt*f_{i} + dt*m_{i}*g_{i}*sum_p(w_{ip})

  where f_{i} is the force in the i face of direction X/Y/Z. Calculated by calcDeviatoricForce
  --------------------------------------------------------------------------------------------------------------
  */

  float BComponent=0.0;

  //Get particle velocity and mass
  float particleMass=_particle->getMass();
  Eigen::Vector3f particleVelocity=_particle->getPreviousVelocity();

  //Add previous velocity contribution
  float velocityComponent=_eVector.dot(particleVelocity);
  BComponent+=(_weight*particleMass*velocityComponent);

  //Add deviatoric force contribution
  BComponent+=(m_dt*_deviatoricForce);

  //Add external force contribution
  float externalForceComponent=_eVector.dot(m_externalForce);
  BComponent+=(m_dt*_massFace*_weight*externalForceComponent);

  //Return result
  return BComponent;

}

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcAComponent_DeviatoricVelocity_New(Particle *_particle, Eigen::Vector3f _weight_i_diff, Eigen::Vector3f _weight_j_diff, Eigen::Vector3f _eVector)
{

}

//----------------------------------------------------------------------------------------------------------------------
