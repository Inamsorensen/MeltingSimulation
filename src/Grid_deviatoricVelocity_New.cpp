#include "Grid.h"

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcDeviatoricContributions(Particle *_particle, int _cellIndex, int _iIndex, int _jIndex, int _kIndex,
                                       float _weightX, float _weightY, float _weightZ)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Calc differentiated weight

  Add to force

  Add to B component

  Loop over all possible cells for each cell

    Add to A component
  ------------------------------------------------------------------------------------------------------
  */

  //Set e_{a(i)} vectors
  Eigen::Vector3f e_x(1.0, 0.0, 0.0);
  Eigen::Vector3f e_y(0.0, 1.0, 0.0);
  Eigen::Vector3f e_z(0.0, 0.0, 1.0);

  //Get particle position
  Eigen::Vector3f particlePosition=_particle->getPosition();

  //Get number of particles in faces
  int noParticles_FaceX=m_cellFacesX[_cellIndex]->m_noParticlesContributing;
  int noParticles_FaceY=m_cellFacesY[_cellIndex]->m_noParticlesContributing;
  int noParticles_FaceZ=m_cellFacesZ[_cellIndex]->m_noParticlesContributing;

  //Get differentiated weights
  Eigen::Vector3f weightDiff_FaceX;
  Eigen::Vector3f weightDiff_FaceY;
  Eigen::Vector3f weightDiff_FaceZ;
  calcWeight_cubicBSpline_Diff(particlePosition, _iIndex, _jIndex, _kIndex, weightDiff_FaceX, weightDiff_FaceY, weightDiff_FaceZ);

  if (noParticles_FaceX>0)
  {
    //Get deviatoric forces
    float force_FaceX=calcDeviatoricForce_New(_particle, e_x, weightDiff_FaceX);

    //Add force to cell faces
    m_cellFacesX[_cellIndex]->m_deviatoricForce+=force_FaceX;

    //Get mass of faces
    float mass_FaceX=m_cellFacesX[_cellIndex]->m_mass;

    //Calculate B Components
    float BComponentX=calcBComponent_DeviatoricVelocity_New(_particle, e_x, _weightX, force_FaceX, mass_FaceX);

    //Add contributions to B vectors
    m_Bvector_deviatoric_X(_cellIndex)+=BComponentX;

  }

  if (noParticles_FaceY>0)
  {
    //Get deviatoric forces
    float force_FaceY=calcDeviatoricForce_New(_particle, e_y, weightDiff_FaceY);

    //Add force to cell faces
    m_cellFacesY[_cellIndex]->m_deviatoricForce+=force_FaceY;

    //Get mass of faces
    float mass_FaceY=m_cellFacesY[_cellIndex]->m_mass;

    //Calculate B Components
    float BComponentY=calcBComponent_DeviatoricVelocity_New(_particle, e_y, _weightY, force_FaceY, mass_FaceY);

    //Add contributions to B vectors
    m_Bvector_deviatoric_Y(_cellIndex)+=BComponentY;

  }

  if (noParticles_FaceZ>0)
  {
    //Get deviatoric forces
    float force_FaceZ=calcDeviatoricForce_New(_particle, e_z, weightDiff_FaceZ);

    //Add force to cell faces
    m_cellFacesZ[_cellIndex]->m_deviatoricForce+=force_FaceZ;

    //Get mass of faces
    float mass_FaceZ=m_cellFacesZ[_cellIndex]->m_mass;

    //Calculate B Components
    float BComponentZ=calcBComponent_DeviatoricVelocity_New(_particle, e_z, _weightZ, force_FaceZ, mass_FaceZ);

    //Add contributions to B vectors
    m_Bvector_deviatoric_Z(_cellIndex)+=BComponentZ;

  }

  if (m_isImplictIntegration==true)
  {
    calcAComponent_DeviatoricVelocity_New(_particle, _cellIndex, weightDiff_FaceX, weightDiff_FaceY, weightDiff_FaceZ);
  }

}

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

void Grid::calcAComponent_DeviatoricVelocity_New(Particle *_particle, int _cellIndex_column, Eigen::Vector3f _weightDiff_FaceX_column, Eigen::Vector3f _weightDiff_FaceY_column, Eigen::Vector3f _weightDiff_FaceZ_column)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------------
  Calculates Ap for for one column in A: the cell at _cellIndex_column

  Determines which cell the particle is in

  Loops over all cells the particle will contribute to

    Calculate weight_diff for each of these

    Multiply calculate the contribution from this particle to the A matrix component at cellIndex_row by cellIndex_column

    Add result to A matrix

  ------------------------------------------------------------------------------------------------------------
  */

  //Set e_{a(i)} vectors
  Eigen::Vector3f e_x(1.0, 0.0, 0.0);
  Eigen::Vector3f e_y(0.0, 1.0, 0.0);
  Eigen::Vector3f e_z(0.0, 0.0, 1.0);

  //Get number of particles in faces
  int noParticles_FaceX_column=m_cellFacesX[_cellIndex_column]->m_noParticlesContributing;
  int noParticles_FaceY_column=m_cellFacesY[_cellIndex_column]->m_noParticlesContributing;
  int noParticles_FaceZ_column=m_cellFacesZ[_cellIndex_column]->m_noParticlesContributing;

  //Get particle parameters necessary: V_{p} and F_{Ep}
  float particleVolume=_particle->getVolume();
  Eigen::Matrix3f deformGradElastic=_particle->getDeformationElastic();
  Eigen::Matrix3f deformGradElastic_trans=deformGradElastic.transpose();

  //Calculate Ap component for the particle and cellIndex_column
  Eigen::Matrix3f ApComponentX;
  Eigen::Matrix3f ApComponentY;
  Eigen::Matrix3f ApComponentZ;

  if (noParticles_FaceX_column>0)
  {
    ApComponentX=calcApComponent_DeviatoricVelocity_New(_particle, _weightDiff_FaceX_column, e_x);
  }

  if (noParticles_FaceY_column>0)
  {
    ApComponentY=calcApComponent_DeviatoricVelocity_New(_particle, _weightDiff_FaceY_column, e_y);
  }

  if (noParticles_FaceZ_column>0)
  {
    ApComponentZ=calcApComponent_DeviatoricVelocity_New(_particle, _weightDiff_FaceZ_column, e_z);
  }

  //Calculate position of particle in grid
  //To calc position of particle, need origin of grid corner, not centre of first grid cell.
  float halfCellSize=m_cellSize/2.0;
  Eigen::Vector3f gridEdgePosition=m_origin;
  gridEdgePosition(0)-=halfCellSize;
  gridEdgePosition(1)-=halfCellSize;
  gridEdgePosition(2)-=halfCellSize;

  //Get particle position in grid
  Eigen::Vector3f particlePosition=_particle->getPosition();
  Eigen::Vector3i particleIndex=MathFunctions::getParticleGridCell(particlePosition, m_cellSize, gridEdgePosition);

  //Loop over i+-2, j+-2, k+-2 to get cells that particle will contribute to
  int iParticle=particleIndex(0);
  int jParticle=particleIndex(1);
  int kParticle=particleIndex(2);

  for (int kIndex_row=(kParticle-2); kIndex_row<(kParticle+4); kIndex_row++)
  {
    for (int jIndex_row=(jParticle-2); jIndex_row<(jParticle+4); jIndex_row++)
    {
      for (int iIndex_row=(iParticle-2); iIndex_row<(iParticle+4); iIndex_row++)
      {
        //Check that not outside grid
        if (iIndex_row>=0 && iIndex_row<=(m_noCells-1) && jIndex_row>=0 && jIndex_row<=(m_noCells-1) && kIndex_row>=0 && kIndex_row<=(m_noCells-1))
        {
          //Get cell index of row
          int cellIndex_row=MathFunctions::getVectorIndex(iIndex_row, jIndex_row, kIndex_row, m_noCells);

          //Get number of particles in faces
          int noParticles_FaceX_row=m_cellFacesX[cellIndex_row]->m_noParticlesContributing;
          int noParticles_FaceY_row=m_cellFacesY[cellIndex_row]->m_noParticlesContributing;
          int noParticles_FaceZ_row=m_cellFacesZ[cellIndex_row]->m_noParticlesContributing;

          //Get differentiated weights
          Eigen::Vector3f weightDiff_FaceX_row;
          Eigen::Vector3f weightDiff_FaceY_row;
          Eigen::Vector3f weightDiff_FaceZ_row;
          calcWeight_cubicBSpline_Diff(particlePosition, iIndex_row, jIndex_row, kIndex_row, weightDiff_FaceX_row, weightDiff_FaceY_row, weightDiff_FaceZ_row);

          //Multiply Ap with row variables
          float AComponentX=0.0;
          float AComponentY=0.0;
          float AComponentZ=0.0;

          if (noParticles_FaceX_column>0 && noParticles_FaceX_row>0)
          {
            Eigen::Vector3f part1_X=deformGradElastic_trans*weightDiff_FaceX_row;
            Eigen::Vector3f part2_X=ApComponentX*part1_X;
            AComponentX=e_x.dot(part2_X);
            AComponentX*=(particleVolume*((pow(m_dt, 2.0))/2.0));
          }

          if (noParticles_FaceY_column>0 && noParticles_FaceY_row>0)
          {
            Eigen::Vector3f part1_Y=deformGradElastic_trans*weightDiff_FaceY_row;
            Eigen::Vector3f part2_Y=ApComponentY*part1_Y;
            AComponentY=e_y.dot(part2_Y);
            AComponentY*=(particleVolume*((pow(m_dt, 2.0))/2.0));
          }

          if (noParticles_FaceZ_column>0 && noParticles_FaceZ_row>0)
          {
            Eigen::Vector3f part1_Z=deformGradElastic_trans*weightDiff_FaceZ_row;
            Eigen::Vector3f part2_Z=ApComponentZ*part1_Z;
            AComponentZ=e_z.dot(part2_Z);
            AComponentZ*=(particleVolume*((pow(m_dt, 2.0))/2.0));
          }


          //Add mass to diagonal elements
          if (cellIndex_row==_cellIndex_column)
          {
            if (noParticles_FaceX_row>0)
            {
              AComponentX+=m_cellFacesX[cellIndex_row]->m_mass;
            }

            if (noParticles_FaceY_row>0)
            {
              AComponentY+=m_cellFacesY[cellIndex_row]->m_mass;
            }

            if (noParticles_FaceZ_row>0)
            {
              AComponentZ+=m_cellFacesZ[cellIndex_row]->m_mass;
            }
          }


          //Add contribution from this particle to components of A matrix
          if (noParticles_FaceX_column>0 && noParticles_FaceX_row>0)
          {
            m_Amatrix_deviatoric_X(cellIndex_row, _cellIndex_column)+=AComponentX;
          }

          if (noParticles_FaceY_column>0 && noParticles_FaceY_row>0)
          {
            m_Amatrix_deviatoric_Y(cellIndex_row, _cellIndex_column)+=AComponentY;
          }

          if (noParticles_FaceZ_column>0 && noParticles_FaceZ_row>0)
          {
            m_Amatrix_deviatoric_Z(cellIndex_row, _cellIndex_column)+=AComponentZ;
          }

        }
      }
    }
  }

}

//----------------------------------------------------------------------------------------------------------------------

Eigen::Matrix3f Grid::calcApComponent_DeviatoricVelocity_New(Particle *_particle, Eigen::Vector3f _weight_diff_column, Eigen::Vector3f _eVector)
{
  /* Outline
  --------------------------------------------------------------------------------------------------------------
    Calculate Z=eVector*weight_diff_trans*deformGradElastic

    Calculate deltaJF=B:Z

    Calculate deltaR

    Calculate C:deltaJF -> d2YdF2:dF

    Calculate C_hat:Z by calculating each part and adding

    Return Ap=C_hat:Z

  --------------------------------------------------------------------------------------------------------------
  */

  //Set up storage for different parts
  Eigen::Matrix3f part1;
  Eigen::Matrix3f part2;
  Eigen::Matrix3f part3;
  Eigen::Matrix3f part4;

  //Get parameters from particle
  Eigen::Matrix3f deformationGradElastic=_particle->getDeformationElastic();
  Eigen::Matrix3f R_deformElastic_Deviatoric=_particle->getR_deformationElastic_Deviatoric();
  Eigen::Matrix3f S_deformElastic_Deviatoric=_particle->getS_deformationElastic_Deviatoric();
  float lameMu=_particle->getLameMu();
  float dimension=_particle->getDimension();
  float detDeformElastic=_particle->getDetDeformationElastic();

  Eigen::Matrix3f deformElastic_trans=deformationGradElastic.transpose();
  Eigen::Matrix3f deformElastic_trans_inverse=deformElastic_trans.inverse();

  //Calculate Z
  Eigen::Matrix3f Z_matrix;
  Z_matrix(0,0)=_eVector(0)*_weight_diff_column(0);
  Z_matrix(0,1)=_eVector(0)*_weight_diff_column(1);
  Z_matrix(0,2)=_eVector(0)*_weight_diff_column(2);
  Z_matrix(1,0)=_eVector(1)*_weight_diff_column(0);
  Z_matrix(1,1)=_eVector(1)*_weight_diff_column(1);
  Z_matrix(1,2)=_eVector(1)*_weight_diff_column(2);
  Z_matrix(2,0)=_eVector(2)*_weight_diff_column(0);
  Z_matrix(2,1)=_eVector(2)*_weight_diff_column(1);
  Z_matrix(2,2)=_eVector(2)*_weight_diff_column(2);

  Z_matrix=(Z_matrix*deformationGradElastic);

  Eigen::Matrix3f Z_matrix_trans=Z_matrix.transpose();


  //Calculate deltaJF=B:Z_matrix
  Eigen::Matrix3f deltaJF=_particle->getDeformEDevDiff_Z(Z_matrix);

  //Calculate first part of differential
  //Calculate deltaR
  Eigen::Matrix3f deltaR=calculate_dR(deltaJF, R_deformElastic_Deviatoric, S_deformElastic_Deviatoric);

  //Calculate d2YdF2
  Eigen::Matrix3f d2YdF2=(2.0*lameMu)*(deltaJF-deltaR);

  //Find part 1 of differential
  part1=_particle->getZ_DeformEDevDiff(d2YdF2);


  //Get dYdF which is used in part 2,3,4
  Eigen::Matrix3f dYdFE=_particle->getPotentialEnergyDiff();

  //Calculate a=-1/dimensions and J^a
  float aConstant=(-1.0/dimension);
  float JaConstant=pow(detDeformElastic, aConstant);

  //Calculate part 2
  Eigen::Matrix3f part2_1=MathFunctions::matrixElementMultiplication(deformElastic_trans_inverse, Z_matrix);
  Eigen::Matrix3f part2_2=_particle->getZ_DeformEDevDiff(dYdFE);
  part2=part2_1*part2_2;
  part2*=aConstant;

  //Calculate part 3
  part3=MathFunctions::matrixElementMultiplication(dYdFE, Z_matrix);
  part3=part3*deformElastic_trans_inverse;
  part3*=(aConstant*JaConstant);

  //Calculate part 4
  part4=MathFunctions::matrixElementMultiplication(dYdFE, deformationGradElastic);
  part4=part4*deformElastic_trans_inverse;
  part4=part4*Z_matrix_trans;
  part4=part4*deformElastic_trans_inverse;
  part4*=(aConstant*JaConstant);
  part4*=(-1.0);


  //Add all parts to find C_hat:Z
  Eigen::Matrix3f Ap_matrix=(part1+part2+part3+part4);

  //Return result
  return Ap_matrix;
}

//----------------------------------------------------------------------------------------------------------------------

Eigen::Matrix3f Grid::calculate_dR_New(const Eigen::Matrix3f &_deltaDeformElastic_Deviatoric, const Eigen::Matrix3f &_R_deformElastic_Deviatoric, const Eigen::Matrix3f &_S_deformElastic_Deviatoric)
{
  /* Outline - For one particle
  ------------------------------------------------------------------------------------------------------------
  Calculate R^T*dJF and dJF^T*R

  Create b vector

  Create A vector from S from polar decomposition

  Set matrix for R^T*dR

  Verify that R^T*dJF-dJF^T*R is anti symmetric. Otherwise the remaining 3 equations are not linearly dependent

  Solve linear system to obtain x, y and z

  Check these against other three functions, or create second matrix and solve for this too and compare

  Put values into R^T*dR

  Multiply by R to obtain dR

  ------------------------------------------------------------------------------------------------------------
  */

  //Calculate R^T*dJF and dJF^T*R
  Eigen::Matrix3f R_transpose=_R_deformElastic_Deviatoric.transpose();
  Eigen::Matrix3f deltaJF_transpose=_deltaDeformElastic_Deviatoric.transpose();

  Eigen::Matrix3f B_matrix=(R_transpose*_deltaDeformElastic_Deviatoric);
  B_matrix-=(deltaJF_transpose*_R_deformElastic_Deviatoric);

  //Create A matrix and b vector
  Eigen::Matrix3f A_matrix;
  Eigen::Vector3f b_vector;
  Eigen::Vector3f solution_vector;

  //Insert elements to A
  ///NB! Double check these elements
  A_matrix(0,0)=_S_deformElastic_Deviatoric(0,0)+_S_deformElastic_Deviatoric(1,1);
  A_matrix(0,1)=_S_deformElastic_Deviatoric(1,2);
  A_matrix(0,2)=(-1.0)*_S_deformElastic_Deviatoric(0,2);
  A_matrix(1,0)=A_matrix(0,1);
  A_matrix(1,1)=_S_deformElastic_Deviatoric(0,0)+_S_deformElastic_Deviatoric(2,2);
  A_matrix(1,2)=_S_deformElastic_Deviatoric(0,1);
  A_matrix(2,0)=A_matrix(0,2);
  A_matrix(2,1)=A_matrix(1,2);
  A_matrix(2,2)=_S_deformElastic_Deviatoric(1,1)+_S_deformElastic_Deviatoric(2,2);

  //Insert elements to b vector
  b_vector(0)=B_matrix(0,1);
  b_vector(1)=B_matrix(0,2);
  b_vector(2)=B_matrix(1,2);

  //Initialise solution to zero
  solution_vector.setZero();

  //Verify that B matrix is antisymmetric. Otherwise non-dependent equations.
  //This check might not work due to floating point inaccuracies

  float tolerance=pow(10.0,-7.0);

  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      //Calculate B(i,j)+B(j,i). Should be zero cause antisymmetric matrices have opposite negatives
      //and zero on diagonal
      float antisymmetryCheck=B_matrix(i,j)+B_matrix(j,i);

      if (antisymmetryCheck>tolerance)
      {
        throw std::invalid_argument("The vector used to calculate delta R is not anti symmetric.");
      }
    }
  }

  //Use linear solver to find x, y, z for RT*deltaR
  MathFunctions::linearSystemSolve(A_matrix, b_vector, solution_vector);


  //Insert solution to RT*deltaR
  Eigen::Matrix3f Rtrans_deltaR;
  Rtrans_deltaR(0,0)=0.0;
  Rtrans_deltaR(0,1)=solution_vector(0);
  Rtrans_deltaR(0,2)=solution_vector(1);
  Rtrans_deltaR(1,0)=(-1.0)*Rtrans_deltaR(0,1);
  Rtrans_deltaR(1,1)=0.0;
  Rtrans_deltaR(1,2)=solution_vector(2);
  Rtrans_deltaR(2,0)=(-1.0)*Rtrans_deltaR(0,2);
  Rtrans_deltaR(2,1)=(-1.0)*Rtrans_deltaR(1,2);
  Rtrans_deltaR(2,2)=0.0;

  //Multiply with R to get deltaR
  Eigen::Matrix3f deltaR=_R_deformElastic_Deviatoric*Rtrans_deltaR;

  //Return result
  return deltaR;
}

//----------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------

void Grid::explicitUpdate_DeviatoricVelocity_New()
{
#pragma omp parallel for
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    m_cellFacesX[cellIndex]->m_velocity=m_Bvector_deviatoric_X(cellIndex);
    m_cellFacesY[cellIndex]->m_velocity=m_Bvector_deviatoric_Y(cellIndex);
    m_cellFacesZ[cellIndex]->m_velocity=m_Bvector_deviatoric_Z(cellIndex);
  }
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::implicitUpdate_DeviatoricVelocity_New()
{

  //Call MINRES to calculate new velocity values
  Eigen::VectorXf solution_X(m_totNoCells);
  solution_X.setZero();
  Eigen::VectorXf solution_Y(m_totNoCells);
  solution_Y.setZero();
  Eigen::VectorXf solution_Z(m_totNoCells);
  solution_Z.setZero();
  Eigen::MatrixXf emptyPreconditioner;

  float shift=(0.0);
  float tolerance=0.0000001;
  int maxNoLoops=20;

//  //Testing symmetry of A matrix
//  Eigen::MatrixXf A_X_trans=m_Amatrix_deviatoric_X.transpose();
//  Eigen::MatrixXf test=m_Amatrix_deviatoric_X-A_X_trans;

//  for (int testItr=0; testItr<m_totNoCells; testItr++)
//  {
//    for (int testItr2=0; testItr2<m_totNoCells; testItr2++)
//    {
//      if (testItr==testItr2)
//      {
//        if ((m_Amatrix_deviatoric_X(testItr, testItr2))>0)
//        {
//          std::cout<<"Index: ["<<testItr<<", "<<testItr2<<"] : "<<m_Amatrix_deviatoric_X(testItr, testItr2)<<"\n";
//        }
//      }
//    }
//  }

//  //Remake matrix to check if determinant is zero
//  Eigen::MatrixXf testSingular;
//  testSingular=testSingular.Identity(m_totNoCells, m_totNoCells);
//  for (int testItr=0; testItr<m_totNoCells; testItr++)
//  {
//    for (int testItr2=0; testItr2<m_totNoCells; testItr2++)
//    {
//      if (testItr!=testItr2)
//      {
//        if (m_Amatrix_deviatoric_X(testItr, testItr2)!=0)
//        {
//          testSingular(testItr,testItr2)=m_Amatrix_deviatoric_X(testItr, testItr2);
//        }
//      }
//    }
//  }
//  float determinant=testSingular.determinant();

  //Set to display details or not
  bool displayDetails=true;
//  bool displayDetails=false;

  //Solve system using MINRES
  MathFunctions::MinRes(m_Amatrix_deviatoric_X, m_Bvector_deviatoric_X, solution_X, emptyPreconditioner, shift, maxNoLoops, tolerance, displayDetails);
  MathFunctions::MinRes(m_Amatrix_deviatoric_Y, m_Bvector_deviatoric_Y, solution_Y, emptyPreconditioner, shift, maxNoLoops, tolerance, displayDetails);
  MathFunctions::MinRes(m_Amatrix_deviatoric_Z, m_Bvector_deviatoric_Z, solution_Z, emptyPreconditioner, shift, maxNoLoops, tolerance, displayDetails);


  //Read in solutions
  #pragma omp parallel for
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    m_cellFacesX[cellIndex]->m_velocity=solution_X[cellIndex];
    m_cellFacesY[cellIndex]->m_velocity=solution_Y[cellIndex];
    m_cellFacesZ[cellIndex]->m_velocity=solution_Z[cellIndex];
  }

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::setBoundaryVelocity()
{
  /* Outline
  ---------------------------------------------------------------------------------------------------------------------
  Loop over all cells
    Check if faces are actually boundaries to colliding cell
    Done by checking whether current cell is colliding. If yes, then all faces are collision boundaries
    If current cell not collidining but faces are, check cells before them
      If cell before doesn't exist cause current cell is ijk=0, then set faces as collision boundaries
      Else check cell before in either face direction and if that cell is colliding, set face to collision boundary

    When face found to be collision boundary, set its velocity to zero for stick.
    This will need to change if object it's colliding with has a velocity.

  ---------------------------------------------------------------------------------------------------------------------
  */

#pragma omp parallel for
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    //Test parallel
//    printf("The parallel region is executed by thread %d\n", omp_get_thread_num());

    //Check whether current cell is colliding, if so set stick collision to all faces
    if (m_cellCentres[cellIndex]->m_state==State::Colliding)
    {
      m_cellFacesX[cellIndex]->m_velocity=0.0;
      m_cellFacesY[cellIndex]->m_velocity=0.0;
      m_cellFacesZ[cellIndex]->m_velocity=0.0;
    }
    else
    {
      //Get cell index in ijk values
      int iIndex=m_cellCentres[cellIndex]->m_iIndex;
      int jIndex=m_cellCentres[cellIndex]->m_jIndex;
      int kIndex=m_cellCentres[cellIndex]->m_kIndex;

      //If current cell isn't colliding, then still need to check the cell before in each direction, unless i||j||k=0
      //FaceX
      if (m_cellFacesX[cellIndex]->m_state==State::Colliding)
      {
        //Set velocity to zero if iIndex=0
        if (iIndex==0)
        {
          m_cellFacesX[cellIndex]->m_velocity=0.0;
        }
        //Else check the cell before it in i direction
        else
        {
          //Get index of cell before in i direction
          int neighbourIndex=MathFunctions::getVectorIndex(iIndex-1, jIndex, kIndex, m_noCells);

          //Check if colliding
          if (m_cellCentres[neighbourIndex]->m_state==State::Colliding)
          {
            m_cellFacesX[cellIndex]->m_velocity=0.0;
          }
        }
      }

      //FaceY
      if (m_cellFacesY[cellIndex]->m_state==State::Colliding)
      {
        //Set velocity to zero if jIndex=0
        if (jIndex==0)
        {
          m_cellFacesY[cellIndex]->m_velocity=0.0;
        }
        //Else check the cell before it in j direction
        else
        {
          //Get index of cell before in j direction
          int neighbourIndex=MathFunctions::getVectorIndex(iIndex, jIndex-1, kIndex, m_noCells);

          //Check if colliding
          if (m_cellCentres[neighbourIndex]->m_state==State::Colliding)
          {
            m_cellFacesY[cellIndex]->m_velocity=0.0;
          }
        }
      }
      //FaceZ
      if (m_cellFacesZ[cellIndex]->m_state==State::Colliding)
      {
        //Set velocity to zero if kIndex=0
        if (kIndex==0)
        {
          m_cellFacesZ[cellIndex]->m_velocity=0.0;
        }
        //Else check the cell before it in k direction
        else
        {
          //Get index of cell before in k direction
          int neighbourIndex=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex-1, m_noCells);

          //Check if colliding
          if (m_cellCentres[neighbourIndex]->m_state==State::Colliding)
          {
            m_cellFacesZ[cellIndex]->m_velocity=0.0;
          }
        }
      }
    }
  }

}

//----------------------------------------------------------------------------------------------------------------------
