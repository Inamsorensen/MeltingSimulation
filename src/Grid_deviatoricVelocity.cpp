#include "Grid.h"

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcDeviatoricForce(Particle* _particle, Eigen::Vector3f _eVector, Eigen::Vector3f _weightDiff)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  f_{i}=-sum_p(V_{p}*e_{a(i)}^T*dYdFE_{p}:B_{kmij}*FE_{p}^T*cubicBSpline_Diff_{ip}

  dYdFE_{p}=2*lame_mu_{p}*JE^a*FE_{p} - 2*lame_mu_{p}*RE_{p} where a=-1/d where d=3? for dimensions and RE is the polar decomposition of JE^a*FE

  B_{kmij}=JE^a*I + a*JE^a*FE_{p}^-T*FE_{p} - Stored in particle

  Loop over all grid cells
    Loop over all particles in grid cell

      calc dYdFE_{p}

      special multiply dYdFE_{p} and deformationElastic_Deviatoric_Diff

      multiply with FE^T

      multiply with V_{p}

      Loop over faces
         Multiply everything with weight and e specific to face

         Add to face's force

  ----------------------------------------------------------------------------------------------------------------
  */

  //Set final result
  float finalResult;

  //Get deformation gradient elastic, FE, and transpose it
  Eigen::Matrix3f deformationElastic=_particle->getDeformationElastic();
  Eigen::Matrix3f deformationElastic_trans=deformationElastic.transpose();

  //Get particle volume
  float particleVolume=_particle->getVolume();

  //Get deformation gradient deviatoric, J^{-1/d}FE
  Eigen::Matrix3f deformationElastic_Deviatoric=_particle->getDeformationElastic_Deviatoric();

  //Get R from polar decomposition of deformation gradient deviatoric, J^{-1/d}FE
  Eigen::Matrix3f R_deformationElastic_Deviatoric=_particle->getR_deformationElastic_Deviatoric();

  //Get Lame Mu coefficient
  float lameMu=_particle->getLameMu();

  //Calculate dYdFE=2*mu*(FE-RE)
  Eigen::Matrix3f dYdFE=(2.0*lameMu)*(deformationElastic_Deviatoric-R_deformationElastic_Deviatoric);

  //Pass in dYdFE to particle's getZ_DeformEDevDiff. Return is dY^{hat}dFE
  Eigen::Matrix3f dY_hat_dFE=_particle->getZ_DeformEDevDiff(dYdFE);

  //Multiply FE^{T} with differentiated weight
  Eigen::Vector3f result_vec1=deformationElastic_trans*_weightDiff;

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

void Grid::calcDeviatoricVelocity()
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Need some of the same variables for calculating deviatoric force and elements of A for implicit velocity calc
  Hence this function loops over all particles inside each grid cell and all grid cells. Other functions accessed
  on a per particle basis

  Loop over grid cells

    Need to set up A and b storage

    Loop over particles in grid cell

      Run calc deviatoric force

      Calc sum weights for external forces

    End loop particles

    Run calc B element

    Input into B vector - one for each face

  End loop grid cells

  Set up A matrix if implicit calc on

    Loop over grid cells
      Loop over grid cells
        Run calc A element
        Input into A vector

  Run implicit solve - one thread for each face?

  Store results
    Loop over grid cells again to store


  TODO: Set boundary velocities to go into b

  ----------------------------------------------------------------------------------------------------------------
  */  

  //Set e_{a(i)} vectors
  Eigen::Vector3f e_x(1.0, 0.0, 0.0);
  Eigen::Vector3f e_y(0.0, 1.0, 0.0);
  Eigen::Vector3f e_z(0.0, 0.0, 1.0);

  //Calc total number of cells
  int totNoCells=pow(m_noCells,3);

  //Set up A and b storage
  Eigen::VectorXf b_X(totNoCells);
  Eigen::VectorXf b_Y(totNoCells);
  Eigen::VectorXf b_Z(totNoCells);

#pragma omp parallel for
  for (int cellIndex=0; cellIndex<totNoCells; cellIndex++)
  {
    //Test parallel
//    printf("The parallel region is executed by thread %d\n", omp_get_thread_num());


    //Store velocities as step n before update them
    m_cellFacesX[cellIndex]->m_previousVelocity=m_cellFacesX[cellIndex]->m_velocity;
    m_cellFacesY[cellIndex]->m_previousVelocity=m_cellFacesY[cellIndex]->m_velocity;
    m_cellFacesZ[cellIndex]->m_previousVelocity=m_cellFacesZ[cellIndex]->m_velocity;


    //Loop over particles in each face of grid cell
    //Calculate force addition and weight addition

    //Face X

    //Set storage variables
    float forceSumX=0.0;
    float weightSumX=0.0;

    int noParticles_FaceX=m_cellFacesX[cellIndex]->m_interpolationData.size();
    for (int particleIterator=0; particleIterator<noParticles_FaceX; particleIterator++)
    {
      //Get cubic B spline
      float weight=m_cellFacesX[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

      //Get cubic B spline differentiated
      Eigen::Vector3f weight_diff=m_cellFacesX[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline_Diff;

      //Get particle pointer
      Particle* particle=m_cellFacesX[cellIndex]->m_interpolationData[particleIterator]->m_particle;

      //Calculate deviatoric force
      float forceFromParticle=calcDeviatoricForce(particle, e_x, weight_diff);

      //Add force to sum
      forceSumX+=forceFromParticle;

      //Add interpolation weight
      weightSumX+=weight;
    }

    //Face Y
    //Set storage variables
    float forceSumY=0.0;
    float weightSumY=0.0;

    int noParticles_FaceY=m_cellFacesY[cellIndex]->m_interpolationData.size();
    for (int particleIterator=0; particleIterator<noParticles_FaceY; particleIterator++)
    {
      //Get cubic B spline
      float weight=m_cellFacesY[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

      //Get cubic B spline differentiated
      Eigen::Vector3f weight_diff=m_cellFacesY[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline_Diff;

      //Get particle pointer
      Particle* particle=m_cellFacesY[cellIndex]->m_interpolationData[particleIterator]->m_particle;

      //Calculate deviatoric force
      float forceFromParticle=calcDeviatoricForce(particle, e_y, weight_diff);

      //Add force to sum
      forceSumY+=forceFromParticle;

      //Add interpolation weight
      weightSumY+=weight;
    }

    //Face Z
    //Set storage variables
    float forceSumZ=0.0;
    float weightSumZ=0.0;

    int noParticles_FaceZ=m_cellFacesZ[cellIndex]->m_interpolationData.size();
    for (int particleIterator=0; particleIterator<noParticles_FaceZ; particleIterator++)
    {
      //Get cubic B spline
      float weight=m_cellFacesZ[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline;

      //Get cubic B spline differentiated
      Eigen::Vector3f weight_diff=m_cellFacesZ[cellIndex]->m_interpolationData[particleIterator]->m_cubicBSpline_Diff;

      //Get particle pointer
      Particle* particle=m_cellFacesZ[cellIndex]->m_interpolationData[particleIterator]->m_particle;

      //Calculate deviatoric force
      float forceFromParticle=calcDeviatoricForce(particle, e_z, weight_diff);

      //Add force to sum
      forceSumZ+=forceFromParticle;

      //Add interpolation weight
      weightSumZ+=weight;
    }

    //Set b components to zero. Will be zero if face has no particles affecting it
    float bComponentX=0.0;
    float bComponentY=0.0;
    float bComponentZ=0.0;

    //Verify that cell face isn't empty, otherwise divide by zero mass
    if (noParticles_FaceX!=0)
    {
      //Calculate b components for each of the faces
      float velocityX=m_cellFacesX[cellIndex]->m_velocity;
      float massX=m_cellFacesX[cellIndex]->m_mass;

      bComponentX=calcBComponent_DeviatoricVelocity(velocityX, massX, forceSumX, weightSumX, e_x);
    }

    if (noParticles_FaceY!=0)
    {
      float velocityY=m_cellFacesY[cellIndex]->m_velocity;
      float massY=m_cellFacesY[cellIndex]->m_mass;

      bComponentY=calcBComponent_DeviatoricVelocity(velocityY, massY, forceSumY, weightSumY, e_y);
    }

    if (noParticles_FaceZ!=0)
    {
      float velocityZ=m_cellFacesZ[cellIndex]->m_velocity;
      float massZ=m_cellFacesZ[cellIndex]->m_mass;

      bComponentZ=calcBComponent_DeviatoricVelocity(velocityZ, massZ, forceSumZ, weightSumZ, e_z);
    }

    //Insert into b vector
    b_X(cellIndex)=bComponentX;
    b_Y(cellIndex)=bComponentY;
    b_Z(cellIndex)=bComponentZ;

    //For now, explicitly update velocity. JUST FOR TESTING
    explicitUpdateVelocity(cellIndex, bComponentX, bComponentY, bComponentZ);

  }

//  implicitUpdateVelocity(b_X, b_Y, b_Z);


}

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcBComponent_DeviatoricVelocity(float _velocity, float _mass, float _deviatoricForce, float _sumWeight, Eigen::Vector3f _eVector)
{
  /* Outline
  --------------------------------------------------------------------------------------------------------------
  b_{i}=v_{i}^{n} + (dt/m_{i})*f_{i} + dt*g_{i}*sum_p(w_{ip})
  --------------------------------------------------------------------------------------------------------------
  */

  //Test parallel
//  printf("B component is calculated by thread %d\n", omp_get_thread_num());

  //Set result to return
  float result;

  //External forces component
  float externalForce=m_externalForce.dot(_eVector);
  externalForce*=(m_dt*_sumWeight);

  //Deviatoric force component
  float deviatoricForce=(m_dt/_mass)*_deviatoricForce;

  //Calculate b component
  result=_velocity+deviatoricForce+externalForce;

  //Return result
  return result;
}

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcAComponent_DeviatoricVelocity(Particle *_particle, Eigen::Vector3f _weight_i_diff, Eigen::Vector3f _weight_j_diff, Eigen::Vector3f _eVector, float _mass_i)
{
  /* Outline
  --------------------------------------------------------------------------------------------------------------
  For each face

    calculate Z

    calculate deltaJF=B:Z

    calculate_dR();

    calculate d2YdFE2

    calculate C_hat:Z by calculating each part and adding

    Multiply by FE^T and weight_diff

    Multiply by e_a^T and V_{p}

    Add to A element for face

  --------------------------------------------------------------------------------------------------------------
  */

  //Set up storage for different parts
  Eigen::Matrix3f part1;
  Eigen::Matrix3f part2;
  Eigen::Matrix3f part3;
  Eigen::Matrix3f part4;

  //Get parameters from particle
  Eigen::Matrix3f deformationGradElastic=_particle->getDeformationElastic();
  Eigen::Matrix3f deformationElastic_Deviatoric=_particle->getDeformationElastic_Deviatoric();
  Eigen::Matrix3f R_deformElastic_Deviatoric=_particle->getR_deformationElastic_Deviatoric();
  Eigen::Matrix3f S_deformElastic_Deviatoric=_particle->getS_deformationElastic_Deviatoric();
  float lameMu=_particle->getLameMu();
  float dimension=_particle->getDimension();
  float detDeformElastic=_particle->getDetDeformationElastic();

  Eigen::Matrix3f deformElastic_trans=deformationGradElastic.transpose();
  Eigen::Matrix3f deformElastic_trans_inverse=deformElastic_trans.inverse();

  //Calculate Z
  Eigen::Matrix3f Z_matrix;
  Z_matrix(0,0)=_eVector(0)*_weight_j_diff(0);
  Z_matrix(0,1)=_eVector(0)*_weight_j_diff(1);
  Z_matrix(0,2)=_eVector(0)*_weight_j_diff(2);
  Z_matrix(1,0)=_eVector(1)*_weight_j_diff(0);
  Z_matrix(1,1)=_eVector(1)*_weight_j_diff(1);
  Z_matrix(1,2)=_eVector(1)*_weight_j_diff(2);
  Z_matrix(2,0)=_eVector(2)*_weight_j_diff(0);
  Z_matrix(2,1)=_eVector(2)*_weight_j_diff(1);
  Z_matrix(2,2)=_eVector(2)*_weight_j_diff(2);

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


  //Calculate dYdF which is used in part 2,3,4
  Eigen::Matrix3f dYdFE=(2.0*lameMu)*(deformationElastic_Deviatoric-R_deformElastic_Deviatoric);


  //Calculate part 2
  Eigen::Matrix3f part2_1=MathFunctions::matrixElementMultiplication(deformElastic_trans_inverse, Z_matrix);
  Eigen::Matrix3f part2_2=_particle->getZ_DeformEDevDiff(dYdFE);
  part2=part2_1*part2_2;
  part2*=(-1.0/dimension);

  //Calculate part 3
  part3=MathFunctions::matrixElementMultiplication(dYdFE, Z_matrix);
  part3=part3*deformElastic_trans_inverse;
  float Ja=pow(detDeformElastic, (-1.0/dimension));
  part3*=((-1.0/dimension)*Ja);

  //Calculate part 4
  part4=MathFunctions::matrixElementMultiplication(dYdFE, deformationGradElastic);
  part4=part4*deformElastic_trans_inverse;
  part4=part4*Z_matrix_trans;
  part4=part4*deformElastic_trans_inverse;
  part4*=((-1.0/dimension)*Ja);
  part4*=(-1.0);


  //Add all parts to find C_hat:Z
  Eigen::Matrix3f Ap_matrix=(part1+part2+part3+part4);


  //Multiply Ap matrix with other particle dependent variables
  Eigen::Matrix3f Ap_deformElastTrans=Ap_matrix*deformElastic_trans;

  Eigen::Vector3f Ap_deformElastTrans_weightDiff=Ap_deformElastTrans*_weight_i_diff;

  float Acomponent=_eVector.dot(Ap_deformElastTrans_weightDiff);

  Acomponent*=_particle->getVolume();

  //Multiply by deltaT^2/2*mass_i
  Acomponent*=((pow(m_dt,2.0))/(2.0*_mass_i));

  //Return result
  return Acomponent;


}

//----------------------------------------------------------------------------------------------------------------------

Eigen::Matrix3f Grid::calculate_dR(const Eigen::Matrix3f &_deltaDeformElastic_Deviatoric, const Eigen::Matrix3f &_R_deformElastic_Deviatoric, const Eigen::Matrix3f &_S_deformElastic_Deviatoric)
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

void Grid::explicitUpdateVelocity(int _cellIndex, float _velocityX, float _velocityY, float _velocityZ)
{
  m_cellFacesX[_cellIndex]->m_velocity=_velocityX;
  m_cellFacesY[_cellIndex]->m_velocity=_velocityY;
  m_cellFacesZ[_cellIndex]->m_velocity=_velocityZ;
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::implicitUpdateVelocity(const Eigen::VectorXf &_bVector_X, const Eigen::VectorXf &_bVector_Y, const Eigen::VectorXf &_bVector_Z)
{
  //Set up face normals
  Eigen::Vector3f e_x(1.0, 0.0, 0.0);
  Eigen::Vector3f e_y(0.0, 1.0, 0.0);
  Eigen::Vector3f e_z(0.0, 0.0, 1.0);

  //Set up A storage
  int totNoCells=pow(m_noCells,3);
  Eigen::MatrixXf A_X(totNoCells, totNoCells);
  Eigen::MatrixXf A_Y(totNoCells, totNoCells);
  Eigen::MatrixXf A_Z(totNoCells, totNoCells);
  A_X.setZero();
  A_Y.setZero();
  A_Z.setZero();

//  printf("Nested parallelism is %s\n", omp_get_nested() ? "supported" : "not supported");

  omp_set_nested(1);
#pragma omp parallel for
  for (int cellIndex_i=0; cellIndex_i<totNoCells; cellIndex_i++)
  {
    //Test parallel
//    printf("Thread %d executes outer parallel region\n", omp_get_thread_num());

    //Calculate number of particles in faces of cellIndex_j
    int noParticles_FaceX_i=m_cellFacesX[cellIndex_i]->m_interpolationData.size();
    int noParticles_FaceY_i=m_cellFacesY[cellIndex_i]->m_interpolationData.size();
    int noParticles_FaceZ_i=m_cellFacesZ[cellIndex_i]->m_interpolationData.size();

    //Loop over cells again to calculate A components
#pragma omp parallel for
    for (int cellIndex_j=0; cellIndex_j<totNoCells; cellIndex_j++)
    {
      //Test parallel
//      printf("Thread %d executes inner parallel region\n", omp_get_thread_num());

      //Calculate number of particles in faces of cellIndex_j
      int noParticles_FaceX_j=m_cellFacesX[cellIndex_j]->m_interpolationData.size();
      int noParticles_FaceY_j=m_cellFacesY[cellIndex_j]->m_interpolationData.size();
      int noParticles_FaceZ_j=m_cellFacesZ[cellIndex_j]->m_interpolationData.size();

      //Face X
      //Loop over particles in cell face X for cellIndex
      for (int particleIterator_i=0; particleIterator_i<noParticles_FaceX_i; particleIterator_i++)
      {
        for (int particleIterator_j=0; particleIterator_j<noParticles_FaceX_j; particleIterator_j++)
        {
          Particle* particleX_i=m_cellFacesX[cellIndex_i]->m_interpolationData[particleIterator_i]->m_particle;
          Particle* particleX_j=m_cellFacesX[cellIndex_j]->m_interpolationData[particleIterator_j]->m_particle;

          if (particleX_i==particleX_j)
          {
            //Calculate Aij component from p

            //Get interpolation weights
            Eigen::Vector3f cubicBSplineDiff_i=m_cellFacesX[cellIndex_i]->m_interpolationData[particleIterator_i]->m_cubicBSpline_Diff;
            Eigen::Vector3f cubicBSplineDiff_j=m_cellFacesX[cellIndex_j]->m_interpolationData[particleIterator_j]->m_cubicBSpline_Diff;

            float mass_faceX_i=m_cellFacesX[cellIndex_i]->m_mass;
            float Acomponent=calcAComponent_DeviatoricVelocity(particleX_i, cubicBSplineDiff_i, cubicBSplineDiff_j, e_x, mass_faceX_i);

            A_X(cellIndex_i, cellIndex_j)+=Acomponent;
          }
        }
      }

      //Face Y
      //Loop over particles in cell face Y for cellIndex
      for (int particleIterator_i=0; particleIterator_i<noParticles_FaceY_i; particleIterator_i++)
      {
        for (int particleIterator_j=0; particleIterator_j<noParticles_FaceY_j; particleIterator_j++)
        {
          Particle* particleY_i=m_cellFacesY[cellIndex_i]->m_interpolationData[particleIterator_i]->m_particle;
          Particle* particleY_j=m_cellFacesY[cellIndex_j]->m_interpolationData[particleIterator_j]->m_particle;

          if (particleY_i==particleY_j)
          {
            //Calculate Aij component from p
            //Get interpolation weights
            Eigen::Vector3f cubicBSplineDiff_i=m_cellFacesY[cellIndex_i]->m_interpolationData[particleIterator_i]->m_cubicBSpline_Diff;
            Eigen::Vector3f cubicBSplineDiff_j=m_cellFacesY[cellIndex_j]->m_interpolationData[particleIterator_j]->m_cubicBSpline_Diff;

            float mass_faceY_i=m_cellFacesY[cellIndex_i]->m_mass;
            float Acomponent=calcAComponent_DeviatoricVelocity(particleY_i, cubicBSplineDiff_i, cubicBSplineDiff_j, e_y, mass_faceY_i);

            A_Y(cellIndex_i, cellIndex_j)+=Acomponent;
          }
        }
      }

      //Face Z
      //Loop over particles in cell face Y for cellIndex
      for (int particleIterator_i=0; particleIterator_i<noParticles_FaceZ_i; particleIterator_i++)
      {
        for (int particleIterator_j=0; particleIterator_j<noParticles_FaceZ_j; particleIterator_j++)
        {
          Particle* particleZ_i=m_cellFacesZ[cellIndex_i]->m_interpolationData[particleIterator_i]->m_particle;
          Particle* particleZ_j=m_cellFacesZ[cellIndex_j]->m_interpolationData[particleIterator_j]->m_particle;

          if (particleZ_i==particleZ_j)
          {
            //Calculate Aij component from p
            //Get interpolation weights
            Eigen::Vector3f cubicBSplineDiff_i=m_cellFacesZ[cellIndex_i]->m_interpolationData[particleIterator_i]->m_cubicBSpline_Diff;
            Eigen::Vector3f cubicBSplineDiff_j=m_cellFacesZ[cellIndex_j]->m_interpolationData[particleIterator_j]->m_cubicBSpline_Diff;

            float mass_faceZ_i=m_cellFacesZ[cellIndex_i]->m_mass;
            float Acomponent=calcAComponent_DeviatoricVelocity(particleZ_i, cubicBSplineDiff_i, cubicBSplineDiff_j, e_z, mass_faceZ_i);

            A_Z(cellIndex_i, cellIndex_j)+=Acomponent;
          }
        }
      }
    }
  }

  omp_set_nested(0);
  
  //Call MINRES to calculate new velocity values
  Eigen::VectorXf solution_X(totNoCells);
  solution_X.setZero();
  Eigen::VectorXf solution_Y(totNoCells);
  solution_Y.setZero();
  Eigen::VectorXf solution_Z(totNoCells);
  solution_Z.setZero();
  Eigen::MatrixXf emptyPreconditioner;
  float shift=(-1.0);
  float tolerance=0.0000001;
  int maxNoLoops=20;

//  Eigen::MatrixXf A_X_trans=A_X.transpose();
//  Eigen::MatrixXf test=A_X-A_X_trans;

  MathFunctions::MinRes(A_X, _bVector_X, solution_X, emptyPreconditioner, shift, maxNoLoops, tolerance, false);
  MathFunctions::MinRes(A_Y, _bVector_Y, solution_Y, emptyPreconditioner, shift, maxNoLoops, tolerance, false);
  MathFunctions::MinRes(A_Z, _bVector_Z, solution_Z, emptyPreconditioner, shift, maxNoLoops, tolerance, false);


  //Read in solutions
  for (int cellIndex=0; cellIndex<totNoCells; cellIndex++)
  {
    m_cellFacesX[cellIndex]->m_velocity=solution_X[cellIndex];
    m_cellFacesY[cellIndex]->m_velocity=solution_Y[cellIndex];
    m_cellFacesZ[cellIndex]->m_velocity=solution_Z[cellIndex];
  }

}

//----------------------------------------------------------------------------------------------------------------------
