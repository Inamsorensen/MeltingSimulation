#include "Grid.h"

void Grid::updateParticleFromGrid(float _velocityContribAlpha, float _tempContribBeta)
{
  /* Outline
  ---------------------------------------------------------------------------------------------------------------------
  Loop over all cells
    Get quadratic stencil and its differential
    Get velocity and previous velocity
    Get quadratic stencil of centre
    Get temperature and previous temp

    Calc PIC velocity for each face
    Calc FLIP velocity for each face
    Add to particle

    Calc velocity gradient contribution for each face
    Add to particle

    Calc PIC temperature for cell centre
    Calc FLIP temperature for cell centre
    Add to particle
  ---------------------------------------------------------------------------------------------------------------------
  */

  //Set e_{a(i)} vectors
  Eigen::Vector3f e_x(1.0, 0.0, 0.0);
  Eigen::Vector3f e_y(0.0, 1.0, 0.0);
  Eigen::Vector3f e_z(0.0, 0.0, 1.0);

  for (int cellIndex=0; cellIndex<pow(m_noCells,3); ++cellIndex)
  {
    //Get velocity and previous velocity of faces
    float velocity_FaceX=m_cellFacesX[cellIndex]->m_velocity;
    float velocity_FaceY=m_cellFacesY[cellIndex]->m_velocity;
    float velocity_FaceZ=m_cellFacesZ[cellIndex]->m_velocity;
    float prevVelocity_FaceX=m_cellFacesX[cellIndex]->m_previousVelocity;
    float prevVelocity_FaceY=m_cellFacesY[cellIndex]->m_previousVelocity;
    float prevVelocity_FaceZ=m_cellFacesZ[cellIndex]->m_previousVelocity;

    //Get temperature and previous temperature
    float temperature=m_cellCentres[cellIndex]->m_temperature;
    float prevTemperature=m_cellCentres[cellIndex]->m_previousTemperature;

    //Face X

    int noParticles_FaceX=m_cellFacesX[cellIndex]->m_interpolationData.size();
    for (int particleIterator=0; particleIterator<noParticles_FaceX; particleIterator++)
    {
      //Get quadratic stencil and its derivative
      float quadStencil=m_cellFacesX[cellIndex]->m_interpolationData[particleIterator]->m_tightQuadStencil;

      //Check that quad stencil isn't zero
      if (quadStencil!=0)
      {
        Eigen::Vector3f quadStencil_Diff=m_cellFacesX[cellIndex]->m_interpolationData[particleIterator]->m_tightQuadStencil_Diff;

        //PIC velocity
        float velocityPIC=velocity_FaceX*quadStencil;

        //FLIP velocity
        float velocityFLIP=(velocity_FaceX-prevVelocity_FaceX)*quadStencil;

        //Velocity contribution
        float velocityContribution=(_velocityContribAlpha*velocityFLIP)+((1.0-_velocityContribAlpha)*velocityPIC);
        Eigen::Vector3f velContribVector=velocityContribution*e_x;

        //Set up velocity gradient contribution
        Eigen::Matrix3f velGradContribution;
        velGradContribution.setZero();
        velGradContribution(0,0)=velocity_FaceX*quadStencil_Diff(0);
        velGradContribution(0,1)=velocity_FaceX*quadStencil_Diff(1);
        velGradContribution(0,2)=velocity_FaceX*quadStencil_Diff(2);

        //Update particle
        Particle* particle=m_cellFacesX[cellIndex]->m_interpolationData[particleIterator]->m_particle;
        particle->addParticleVelocity(velContribVector);
        particle->addParticleVelocityGradient(velGradContribution);
      }

    }

    //Face Y
    int noParticles_FaceY=m_cellFacesY[cellIndex]->m_interpolationData.size();
    for (int particleIterator=0; particleIterator<noParticles_FaceY; particleIterator++)
    {
      //Get quadratic stencil and its derivative
      float quadStencil=m_cellFacesY[cellIndex]->m_interpolationData[particleIterator]->m_tightQuadStencil;

      if (quadStencil!=0)
      {
        Eigen::Vector3f quadStencil_Diff=m_cellFacesY[cellIndex]->m_interpolationData[particleIterator]->m_tightQuadStencil_Diff;

        //PIC velocity
        float velocityPIC=velocity_FaceY*quadStencil;

        //FLIP velocity
        float velocityFLIP=(velocity_FaceY-prevVelocity_FaceY)*quadStencil;

        //Velocity contribution
        float velocityContribution=(_velocityContribAlpha*velocityFLIP)+((1.0-_velocityContribAlpha)*velocityPIC);
        Eigen::Vector3f velContribVector=velocityContribution*e_y;

        //Set up velocity gradient contribution
        Eigen::Matrix3f velGradContribution;
        velGradContribution.setZero();
        velGradContribution(1,0)=velocity_FaceY*quadStencil_Diff(0);
        velGradContribution(1,1)=velocity_FaceY*quadStencil_Diff(1);
        velGradContribution(1,2)=velocity_FaceY*quadStencil_Diff(2);

        //Update particle
        Particle* particle=m_cellFacesY[cellIndex]->m_interpolationData[particleIterator]->m_particle;
        particle->addParticleVelocity(velContribVector);
        particle->addParticleVelocityGradient(velGradContribution);
      }
    }

    //Face Z
    int noParticles_FaceZ=m_cellFacesZ[cellIndex]->m_interpolationData.size();
    for (int particleIterator=0; particleIterator<noParticles_FaceZ; particleIterator++)
    {
      //Get quadratic stencil and its derivative
      float quadStencil=m_cellFacesZ[cellIndex]->m_interpolationData[particleIterator]->m_tightQuadStencil;

      if (quadStencil!=0)
      {
        Eigen::Vector3f quadStencil_Diff=m_cellFacesZ[cellIndex]->m_interpolationData[particleIterator]->m_tightQuadStencil_Diff;

        //PIC velocity
        float velocityPIC=velocity_FaceZ*quadStencil;

        //FLIP velocity
        float velocityFLIP=(velocity_FaceZ-prevVelocity_FaceZ)*quadStencil;

        //Velocity contribution
        float velocityContribution=(_velocityContribAlpha*velocityFLIP)+((1.0-_velocityContribAlpha)*velocityPIC);
        Eigen::Vector3f velContribVector=velocityContribution*e_z;

        //Set up velocity gradient contribution
        Eigen::Matrix3f velGradContribution;
        velGradContribution.setZero();
        velGradContribution(2,0)=velocity_FaceZ*quadStencil_Diff(0);
        velGradContribution(2,1)=velocity_FaceZ*quadStencil_Diff(1);
        velGradContribution(2,2)=velocity_FaceZ*quadStencil_Diff(2);

        //Update particle
        Particle* particle=m_cellFacesZ[cellIndex]->m_interpolationData[particleIterator]->m_particle;
        particle->addParticleVelocity(velContribVector);
        particle->addParticleVelocityGradient(velGradContribution);
      }
    }

    //Cell centre
    int noParticles_cellCentre=m_cellCentres[cellIndex]->m_interpolationData.size();
    for (int particleIterator=0; particleIterator<noParticles_cellCentre; particleIterator++)
    {
      //Get quadratic stencil
      float quadStencil=m_cellCentres[cellIndex]->m_interpolationData[particleIterator]->m_tightQuadStencil;

      if (quadStencil!=0)
      {
        //PIC temperature
        float temperaturePIC=temperature*quadStencil;

        //FLIP temperature
        float temperatureFLIP=(temperature-prevTemperature)*quadStencil;

        //Calculate temperature contribution
        float temperatureContribution=(_tempContribBeta*temperatureFLIP)+((1.0-_tempContribBeta)*temperaturePIC);

        //Update particle
        Particle* particle=m_cellCentres[cellIndex]->m_interpolationData[particleIterator]->m_particle;
        particle->addParticleTemperature(temperatureContribution);
      }
    }

  }

}


//----------------------------------------------------------------------------------------------------------------------
