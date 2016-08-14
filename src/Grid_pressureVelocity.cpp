#include "Grid.h"

//----------------------------------------------------------------------------------------------------------------------

void Grid::projectVelocity()
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Calculate face densities

  Set up eigen matrices to store values in

  Set up B

  Set up A

  Solve using conjugate gradient

  Set result to cells

  ----------------------------------------------------------------------------------------------------------------
  */

  //Calculate cell face densities for interior cells
//#pragma omp parallel for
//  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
//  {
//    if (m_cellCentres[cellIndex]->m_state==State::Interior)
//    {
//      calcFaceDensities(cellIndex);
//    }
//  }

//Calculate cell face densities for interior cells
#pragma omp parallel for
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    //Only update if cell centres are interior
    if (m_cellCentres[cellIndex]->m_state==State::Interior)
    {
      //Get index of cell
      int iIndex=m_cellCentres[cellIndex]->m_iIndex;
      int jIndex=m_cellCentres[cellIndex]->m_jIndex;
      int kIndex=m_cellCentres[cellIndex]->m_kIndex;

      //Get mass of cell faces
      float massX=m_cellFacesX[cellIndex]->m_mass;
      float massY=m_cellFacesY[cellIndex]->m_mass;
      float massZ=m_cellFacesZ[cellIndex]->m_mass;

      //Calculate control volumes of cell faces
      float volumeX=calcFaceVolume(iIndex, jIndex, kIndex, 0);
      float volumeY=calcFaceVolume(iIndex, jIndex, kIndex, 1);
      float volumeZ=calcFaceVolume(iIndex, jIndex, kIndex, 2);

      //Set cell face densities
      m_cellFacesX[cellIndex]->m_density=massX/volumeX;
      m_cellFacesY[cellIndex]->m_density=massY/volumeY;
      m_cellFacesZ[cellIndex]->m_density=massZ/volumeZ;

      //Check if cells of the upper faces are empty or colliding, if so calculate their density too
      int cellIndex_i1jk=MathFunctions::getVectorIndex(iIndex+1, jIndex, kIndex, m_noCells);
      int cellIndex_ij1k=MathFunctions::getVectorIndex(iIndex, jIndex+1, kIndex, m_noCells);
      int cellIndex_ijk1=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex+1, m_noCells);

      if (m_cellCentres[cellIndex_i1jk]->m_state!=State::Interior)
      {
        float mass=m_cellFacesX[cellIndex_i1jk]->m_mass;
        float volume=calcFaceVolume(iIndex+1, jIndex, kIndex, 0);
        //Make sure volume isn't zero
        if (volume!=0.0)
        {
          m_cellFacesX[cellIndex_i1jk]->m_density=mass/volume;
        }
      }
      if (m_cellCentres[cellIndex_ij1k]->m_state!=State::Interior)
      {
        float mass=m_cellFacesY[cellIndex_ij1k]->m_mass;
        float volume=calcFaceVolume(iIndex, jIndex+1, kIndex, 1);
        //Make sure volume isn't zero
        if (volume!=0.0)
        {
          m_cellFacesY[cellIndex_ij1k]->m_density=mass/volume;
        }
      }
      if (m_cellCentres[cellIndex_ijk1]->m_state!=State::Interior)
      {
        float mass=m_cellFacesZ[cellIndex_ijk1]->m_mass;
        float volume=calcFaceVolume(iIndex, jIndex, kIndex+1, 2);
        //Make sure volume isn't zero
        if (volume!=0.0)
        {
          m_cellFacesZ[cellIndex_ijk1]->m_density=mass/volume;
        }
      }

    }
  }

  //Set up matrices for linear system
  Eigen::SparseMatrix<double> A_matrix(m_totNoCells, m_totNoCells);
  Eigen::VectorXd B_vector(m_totNoCells);
  Eigen::VectorXd solution(m_totNoCells);

  Eigen::MatrixXf A_matrix_test(m_totNoCells, m_totNoCells);
  A_matrix_test.setZero();

  //Initialise A and B to zero
  A_matrix.setZero();
  B_vector.setZero();
  solution.setZero();

  //Calculate A and B elements
  ///NB! This doesn't work when threaded because of the workings of the sparse matrix.
  /// TODO: Check if faster if create just Eigen::MatrixXf then copy values across, or do B_vector separately and threaded
//#pragma omp parallel for
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    //Only fill in interior cells
    if (m_cellCentres[cellIndex]->m_state==State::Interior)
    {
      //Get ijk
      int iIndex=m_cellCentres[cellIndex]->m_iIndex;
      int jIndex=m_cellCentres[cellIndex]->m_jIndex;
      int kIndex=m_cellCentres[cellIndex]->m_kIndex;

      //Calculate B element
      B_vector(cellIndex)=calcBComponent_projectVelocity(cellIndex, iIndex, jIndex, kIndex);

      //Insert A elements
      calcAComponent_projectVelocity(cellIndex, iIndex, jIndex, kIndex, A_matrix, A_matrix_test);
    }

  }


  //Solve system
  float maxLoops=3000;
  float minResidual=0.00001;
  MathFunctions::conjugateGradient(A_matrix, B_vector, solution, maxLoops, minResidual);


  //Use results to calculate projected velocities
#pragma omp parallel for
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    //Only correct faces surrounding interior cells
    if (m_cellCentres[cellIndex]->m_state==State::Interior)
    {
      //Get indices of cell
      int iIndex=m_cellCentres[cellIndex]->m_iIndex;
      int jIndex=m_cellCentres[cellIndex]->m_jIndex;
      int kIndex=m_cellCentres[cellIndex]->m_kIndex;

      //Get indices of surrounding cells
      int indexCell_i1jk=MathFunctions::getVectorIndex(iIndex+1, jIndex, kIndex, m_noCells);
      int indexCell_i_1jk=MathFunctions::getVectorIndex(iIndex-1, jIndex, kIndex, m_noCells);
      int indexCell_ij1k=MathFunctions::getVectorIndex(iIndex, jIndex+1, kIndex, m_noCells);
      int indexCell_ij_1k=MathFunctions::getVectorIndex(iIndex, jIndex-1, kIndex, m_noCells);
      int indexCell_ijk1=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex+1, m_noCells);
      int indexCell_ijk_1=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex-1, m_noCells);

      ///Update all faces or only the non-colliding?

      //Get pressure for all indices
      float pressure_ijk=solution(cellIndex);
//      float pressure_i1jk=solution(indexCell_i1jk);
//      float pressure_i_1jk=solution(indexCell_i_1jk);
//      float pressure_ij1k=solution(indexCell_ij1k);
//      float pressure_ij_1k=solution(indexCell_ij_1k);
//      float pressure_ijk1=solution(indexCell_ijk1);
//      float pressure_ijk_1=solution(indexCell_ijk_1);

        float pressure_i1jk=0.0;
        float pressure_i_1jk=0.0;
        float pressure_ij1k=0.0;
        float pressure_ij_1k=0.0;
        float pressure_ijk1=0.0;
        float pressure_ijk_1=0.0;

      //Get state of neighbouring cells
      State state_i1jk=m_cellCentres[indexCell_i1jk]->m_state;
      State state_i_1jk=m_cellCentres[indexCell_i_1jk]->m_state;
      State state_ij1k=m_cellCentres[indexCell_ij1k]->m_state;
      State state_ij_1k=m_cellCentres[indexCell_ij_1k]->m_state;
      State state_ijk1=m_cellCentres[indexCell_ijk1]->m_state;
      State state_ijk_1=m_cellCentres[indexCell_ijk_1]->m_state;

      //Enforce boundaries
      //If neighbour cells is interior, set to solution
      //if colliding set to pressure of ijk so pressure gradient is zero
      //If empty, leave pressure as zero.
      if (state_i1jk==State::Interior)
      {
        pressure_i1jk=solution(indexCell_i1jk);
      }
      else if (state_i1jk==State::Colliding)
      {
//        pressure_i1jk=pressure_ijk;
        pressure_i1jk=pressure_ijk*MathFunctions::signFunction(pressure_ijk);
      }

      if (state_i_1jk==State::Interior)
      {
        pressure_i_1jk=solution(indexCell_i_1jk);
      }
      else if (state_i_1jk==State::Colliding)
      {
//        pressure_i_1jk=pressure_ijk;
        pressure_i_1jk=pressure_ijk*MathFunctions::signFunction(pressure_ijk);
      }

      if (state_ij1k==State::Interior)
      {
        pressure_ij1k=solution(indexCell_ij1k);
      }
      else if (state_ij1k==State::Colliding)
      {
//        pressure_ij1k=pressure_ijk;
        pressure_ij1k=pressure_ijk*MathFunctions::signFunction(pressure_ijk);
      }

      if (state_ij_1k==State::Interior)
      {
        pressure_ij_1k=solution(indexCell_ij_1k);
      }
      else if (state_ij_1k==State::Colliding)
      {
//        pressure_ij_1k=pressure_ijk;
        pressure_ij_1k=pressure_ijk*MathFunctions::signFunction(pressure_ijk);
      }

      if (state_ijk1==State::Interior)
      {
        pressure_ijk1=solution(indexCell_ijk1);
      }
      else if (state_ijk1==State::Colliding)
      {
//        pressure_ijk1=pressure_ijk;
        pressure_ijk1=pressure_ijk*MathFunctions::signFunction(pressure_ijk);
      }

      if (state_ijk_1==State::Interior)
      {
        pressure_ijk_1=solution(indexCell_ijk_1);
      }
      else if (state_ijk_1==State::Colliding)
      {
//        pressure_ijk_1=pressure_ijk;
        pressure_ijk_1=pressure_ijk*MathFunctions::signFunction(pressure_ijk);
      }


//      //Check pressure isn't zero
//      if (pressure_ijk<0)
//      {
//        pressure_ijk=0.0;
////        pressure_ijk=std::abs(pressure_ijk);
//      }

//      if (pressure_i1jk<0)
//      {
//        pressure_i1jk=0.0;
////        pressure_i1jk=std::abs(pressure_i1jk);
//      }
//      if (pressure_i_1jk<0)
//      {
//        pressure_i_1jk=0.0;
////        pressure_i_1jk=std::abs(pressure_i_1jk);
//      }
//      if (pressure_ij1k<0)
//      {
//        pressure_ij1k=0.0;
////        pressure_ij1k=std::abs(pressure_ij1k);
//      }
//      if (pressure_ij_1k<0)
//      {
//        pressure_ij_1k=0.0;
////        pressure_ij_1k=std::abs(pressure_ij_1k);
//      }
//      if (pressure_ijk1<0)
//      {
//        pressure_ijk1=0.0;
////        pressure_ijk1=std::abs(pressure_ijk1);
//      }
//      if (pressure_ijk_1<0)
//      {
//        pressure_ijk_1=0.0;
////        pressure_ijk_1=std::abs(pressure_ijk_1);
//      }

      //Calculate pressure gradients
      float pressureGradient_i1jk=pressure_i1jk-pressure_ijk;
      float pressureGradient_i_1jk=pressure_ijk-pressure_i_1jk;
      float pressureGradient_ij1k=pressure_ij1k-pressure_ijk;
      float pressureGradient_ij_1k=pressure_ijk-pressure_ij_1k;
      float pressureGradient_ijk1=pressure_ijk1-pressure_ijk;
      float pressureGradient_ijk_1=pressure_ijk-pressure_ijk_1;

      //Calculate constant to be multiplied with gradient
      float constant=m_dt/m_cellSize;
      float constantX=constant/m_cellFacesX[cellIndex]->m_density;
      float constantY=constant/m_cellFacesY[cellIndex]->m_density;
      float constantZ=constant/m_cellFacesZ[cellIndex]->m_density;
//      float constantX=constant;
//      float constantY=constant;
//      float constantZ=constant;


      //Calculate projected velocity
      m_cellFacesX[cellIndex]->m_velocity=m_cellFacesX[cellIndex]->m_velocity - (constantX*pressureGradient_i_1jk);
      m_cellFacesY[cellIndex]->m_velocity=m_cellFacesY[cellIndex]->m_velocity - (constantY*pressureGradient_ij_1k);
      m_cellFacesZ[cellIndex]->m_velocity=m_cellFacesZ[cellIndex]->m_velocity - (constantZ*pressureGradient_ijk_1);

      //Update upper surrounding faces, if the cells they belong to are empty or colliding
      //Face X
      if (m_cellCentres[indexCell_i1jk]->m_state!=State::Interior)
      {
        float constantX1=constant/m_cellFacesX[indexCell_i1jk]->m_density;
        m_cellFacesX[indexCell_i1jk]->m_velocity=m_cellFacesX[indexCell_i1jk]->m_velocity - (constantX1*pressureGradient_i1jk);
      }
      //Face Y
      if (m_cellCentres[indexCell_ij1k]->m_state!=State::Interior)
      {
        float constantY1=constant/m_cellFacesY[indexCell_ij1k]->m_density;
        m_cellFacesY[indexCell_ij1k]->m_velocity=m_cellFacesY[indexCell_ij1k]->m_velocity - (constantY1*pressureGradient_ij1k);
      }
      //Face Z
      if (m_cellCentres[indexCell_ijk1]->m_state!=State::Interior)
      {
        float constantZ1=constant/m_cellFacesZ[indexCell_ijk1]->m_density;
        m_cellFacesZ[indexCell_ijk1]->m_velocity=m_cellFacesZ[indexCell_ijk1]->m_velocity - (constantZ1*pressureGradient_ijk1);
      }

    }

  }


}

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcBComponent_projectVelocity(int _cellIndex, int _iIndex, int _jIndex, int _kIndex)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------

  Calculate constant which will be added to B component:
    constant=-(J_Ec-1)/(dt*J_Ec)

  Get velocities of faces surrounding cell

  Calculate central difference along each direction for cell in question
    central difference X = velocity[i+1/2,j,k] - velocity[i-1/2,j,k]

  Add central differences and divide by cell size to obtain central gradient stencil

  Add constant and return as B component

  ----------------------------------------------------------------------------------------------------------------
  */

  //Get density
  float densityX=m_cellFacesX[_cellIndex]->m_density;
  float densityY=m_cellFacesY[_cellIndex]->m_density;
  float densityZ=m_cellFacesZ[_cellIndex]->m_density;

  float sumDensity=densityX+densityY+densityZ;

    //Calculate constant
    float detDeformationGradElastic=m_cellCentres[_cellIndex]->m_detDeformationGradElastic;
    float constant=(detDeformationGradElastic-1.0);
    constant*=(1.0/(m_dt*detDeformationGradElastic));
    constant*=(-1.0);

    //Add in density here
    constant*=sumDensity;

    ///TO DO: Put in a check to make sure we aren't searching outside grid.
    /// Think this is prevented by only calculating B for interior cells though

    //Get index of i+1,j,k and i,j+1,k and i,j,k+1
    int index_i1jk=MathFunctions::getVectorIndex(_iIndex+1, _jIndex, _kIndex, m_noCells);
    int index_ij1k=MathFunctions::getVectorIndex(_iIndex, _jIndex+1, _kIndex, m_noCells);
    int index_ijk1=MathFunctions::getVectorIndex(_iIndex, _jIndex, _kIndex+1, m_noCells);

    //Get face velocities for all faces surrounding cell centre
    float velocityX_forward=m_cellFacesX[index_i1jk]->m_velocity;
    float velocityX_backward=m_cellFacesX[_cellIndex]->m_velocity;
    float velocityY_forward=m_cellFacesY[index_ij1k]->m_velocity;
    float velocityY_backward=m_cellFacesY[_cellIndex]->m_velocity;
    float velocityZ_forward=m_cellFacesZ[index_ijk1]->m_velocity;
    float velocityZ_backward=m_cellFacesZ[_cellIndex]->m_velocity;

    //Calculate central gradient stencil
//    float centralGradient=(1.0/m_cellSize)*((velocityX_forward-velocityX_backward)+(velocityY_forward-velocityY_backward)+(velocityZ_forward-velocityZ_backward));
    float centralGradient=(1.0/m_cellSize)*((densityX*(velocityX_forward-velocityX_backward))
                                            +(densityY*(velocityY_forward-velocityY_backward))
                                            +(densityZ*(velocityZ_forward-velocityZ_backward)));

    //Calculate B component
    float BComponent=constant-centralGradient;

    return BComponent;

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcAComponent_projectVelocity(int _cellIndex, int _iIndex, int _jIndex, int _kIndex, Eigen::SparseMatrix<double> &o_A, Eigen::MatrixXf &o_A_test)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------

  Get indices of surrounding cells

  Get density of faces belonging to current cell

  Calculate A components
    A[ijk,ijk] = -6.0*(dt/cellSize^2)*((1/densityX) + (1/desityY) + (1/densityZ))

    Unless any of the surrounding cells are colliding, in which case A[ijk,ijk] is reduced by one to enforce
    the Neumann boundary condition

    The surrounding A components are set to
    (dt/cellSize^2)*(1/density_in_direction)    if cell is interior
    0                                           if cell is empty
    0 and subtract from A[ijk,ijk]              if cell is colliding

  Calculate pressure constant
    J_Pc/(J_Ec*lambda*dt)

  Add pressure constant to A[ijk,ijk]

  Insert all A components into matrix

  ----------------------------------------------------------------------------------------------------------------
  */

  //Get indices of surrounding cells
  int cellIndex_i1jk=MathFunctions::getVectorIndex(_iIndex+1, _jIndex, _kIndex, m_noCells);
  int cellIndex_i_1jk=MathFunctions::getVectorIndex(_iIndex-1, _jIndex, _kIndex, m_noCells);
  int cellIndex_ij1k=MathFunctions::getVectorIndex(_iIndex, _jIndex+1, _kIndex, m_noCells);
  int cellIndex_ij_1k=MathFunctions::getVectorIndex(_iIndex, _jIndex-1, _kIndex, m_noCells);
  int cellIndex_ijk1=MathFunctions::getVectorIndex(_iIndex, _jIndex, _kIndex+1, m_noCells);
  int cellIndex_ijk_1=MathFunctions::getVectorIndex(_iIndex, _jIndex, _kIndex-1, m_noCells);

  //Calculate constant=dt/cellSize^2
  float constant=(m_dt/(pow(m_cellSize,2)));

  //Get densities for cell faces of current cell
  float densityX_ijk=m_cellFacesX[_cellIndex]->m_density;
  float densityY_ijk=m_cellFacesY[_cellIndex]->m_density;
  float densityZ_ijk=m_cellFacesZ[_cellIndex]->m_density;

  //Set up sumInvDensity
//  float A_ijk_X=(-2.0);
//  float A_ijk_Y=(-2.0);
//  float A_ijk_Z=(-2.0);
  float A_ijk_X=(2.0);
  float A_ijk_Y=(2.0);
  float A_ijk_Z=(2.0);

  //Initialise A matrix elements
  float A_i1jk=0.0;
  float A_i_1jk=0.0;
  float A_ij1k=0.0;
  float A_ij_1k=0.0;
  float A_ijk1=0.0;
  float A_ijk_1=0.0;
  float A_ijk=0.0;


  //Calculate A element for surrounding cells
  switch (m_cellCentres[cellIndex_i1jk]->m_state)
  {
  case State::Colliding :
  {
//    A_ijk_X+=1.0;
    A_ijk_X-=1.0;
    break;
  }
  case State::Interior :
  {
//    A_i1jk=(constant/densityX_ijk);
//    A_i1jk=constant;
    A_i1jk=(-1.0*constant);
    break;
  }
  default:
  {
    break;
  }
  }

  switch (m_cellCentres[cellIndex_i_1jk]->m_state)
  {
  case State::Colliding :
  {
//    A_ijk_X+=1.0;
    A_ijk_X-=1.0;
    break;
  }
  case State::Interior :
  {
//    A_i_1jk=(constant/densityX_ijk);
//    A_i_1jk=constant;
    A_i_1jk=(-1.0*constant);
    break;
  }
  default:
  {
    break;
  }
  }

  switch (m_cellCentres[cellIndex_ij1k]->m_state)
  {
  case State::Colliding :
  {
//    A_ijk_Y+=1.0;
    A_ijk_Y-=1.0;
    break;
  }
  case State::Interior :
  {
//    A_ij1k=(constant/densityY_ijk);
//    A_ij1k=constant;
    A_ij1k=(-1.0*constant);
    break;
  }
  default:
  {
    break;
  }
  }

  switch (m_cellCentres[cellIndex_ij_1k]->m_state)
  {
  case State::Colliding :
  {
//    A_ijk_Y+=1.0;
    A_ijk_Y-=1.0;
    break;
  }
  case State::Interior :
  {
//    A_ij_1k=(constant/densityY_ijk);
//    A_ij_1k=constant;
    A_ij_1k=(-1.0*constant);
    break;
  }
  default:
  {
    break;
  }
  }

  switch (m_cellCentres[cellIndex_ijk1]->m_state)
  {
  case State::Colliding :
  {
//    A_ijk_Z+=1.0;
    A_ijk_Z-=1.0;
    break;
  }
  case State::Interior :
  {
//    A_ijk1=(constant/densityZ_ijk);
//    A_ijk1=constant;
    A_ijk1=(-1.0*constant);
    break;
  }
  default:
  {
    break;
  }
  }

  switch (m_cellCentres[cellIndex_ijk_1]->m_state)
  {
  case State::Colliding :
  {
//    A_ijk_Z+=1.0;
    A_ijk_Z-=1.0;
    break;
  }
  case State::Interior :
  {
//    A_ijk_1=(constant/densityZ_ijk);
//    A_ijk_1=constant;
    A_ijk_1=(-1.0*constant);
    break;
  }
  default:
  {
    break;
  }
  }

  //Calculate A_ijk from X,Y,Z contributions, density and constant
//  A_ijk=((A_ijk_X/densityX_ijk)+(A_ijk_Y/densityY_ijk)+(A_ijk_Z/densityZ_ijk));
  A_ijk=(A_ijk_X + A_ijk_Y + A_ijk_Z);
  A_ijk*=constant;

  //Add pressureConstant to A_ijk
  float pressureConst;
  float detDeformGradPlastic=m_cellCentres[_cellIndex]->m_detDeformationGradPlastic;
  float detDeformGradElastic=m_cellCentres[_cellIndex]->m_detDeformationGradElastic;
  float lambdaInv=m_cellCentres[_cellIndex]->m_lameLambdaInverse;
  pressureConst=detDeformGradPlastic/detDeformGradElastic;
  pressureConst*=lambdaInv;
  pressureConst*=(1.0/m_dt);

  //Add density here
  pressureConst*=(densityX_ijk + densityY_ijk + densityZ_ijk);

  A_ijk+=pressureConst;


  //Insert values into matrix
  //Might need to do this in main function for parallelisation
  o_A.insert(_cellIndex, _cellIndex)=A_ijk;
  o_A.insert(_cellIndex, cellIndex_i1jk)=A_i1jk;
  o_A.insert(_cellIndex, cellIndex_i_1jk)=A_i_1jk;
  o_A.insert(_cellIndex, cellIndex_ij1k)=A_ij1k;
  o_A.insert(_cellIndex, cellIndex_ij_1k)=A_ij_1k;
  o_A.insert(_cellIndex, cellIndex_ijk1)=A_ijk1;
  o_A.insert(_cellIndex, cellIndex_ijk_1)=A_ijk_1;

  //Insert values into matrix
  //TEST
  o_A_test(_cellIndex, _cellIndex)=A_ijk;
  o_A_test(_cellIndex, cellIndex_i1jk)=A_i1jk;
  o_A_test(_cellIndex, cellIndex_i_1jk)=A_i_1jk;
  o_A_test(_cellIndex, cellIndex_ij1k)=A_ij1k;
  o_A_test(_cellIndex, cellIndex_ij_1k)=A_ij_1k;
  o_A_test(_cellIndex, cellIndex_ijk1)=A_ijk1;
  o_A_test(_cellIndex, cellIndex_ijk_1)=A_ijk_1;


}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcFaceDensities(int _cellIndex)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Leave for now I think

  ----------------------------------------------------------------------------------------------------------------
  */

  //Get cell indices
  int iIndex=m_cellCentres[_cellIndex]->m_iIndex;
  int jIndex=m_cellCentres[_cellIndex]->m_jIndex;
  int kIndex=m_cellCentres[_cellIndex]->m_kIndex;

  //Set storage for volume
  float xFaceVolume=0.0;
  float yFaceVolume=0.0;
  float zFaceVolume=0.0;

  //Set min and max increments for ijk
  int iMinIncrement=-2;
  int jMinIncrement=-2;
  int kMinIncrement=-2;
  int iMaxIncrement=2;
  int jMaxIncrement=2;
  int kMaxIncrement=2;

  //Make sure doesn't loop over outside cells.
  //Sort of makes the non-existing cells seem like collision cells as will add nothing to volume
  // i direction
  if (iIndex==0)
  {
    iMinIncrement=0;
  }
  if (iIndex==1)
  {
    iMinIncrement=-1;
  }
  if (iIndex==(m_noCells-2))
  {
    iMaxIncrement=1;
  }
  if (iIndex==(m_noCells-1))
  {
    iMaxIncrement=0;
  }

  //j direction
  if (jIndex==0)
  {
    jMinIncrement=0;
  }
  if (jIndex==1)
  {
    jMinIncrement=-1;
  }
  if (jIndex==(m_noCells-2))
  {
    jMaxIncrement=1;
  }
  if (jIndex==(m_noCells-1))
  {
    jMaxIncrement=0;
  }

  // k direction
  if (kIndex==0)
  {
    kMinIncrement=0;
  }
  if (kIndex==1)
  {
    kMinIncrement=-1;
  }
  if (kIndex==(m_noCells-2))
  {
    kMaxIncrement=1;
  }
  if (kIndex==(m_noCells-1))
  {
    kMaxIncrement=0;
  }


  //Loop over neighbouring cells that will give cubicBSpline_Integ!=0
  for (int kIndexIncrement=kMinIncrement; kIndexIncrement<(kMaxIncrement+1); kIndexIncrement++)
  {
    for (int jIndexIncrement=jMinIncrement; jIndexIncrement<(jMaxIncrement+1); jIndexIncrement++)
    {
      for (int iIndexIncrement=iMinIncrement; iIndexIncrement<(iMaxIncrement+1); iIndexIncrement++)
      {

        //Get index of incremented cell
        int cellIndexIncremented=MathFunctions::getVectorIndex((iIndex+iIndexIncrement), (jIndex+jIndexIncrement), (kIndex+kIndexIncrement), m_noCells);

        //Check if cell is colliding. If colliding then will add nothing to cell face volume
        if (m_cellCentres[cellIndexIncremented]->m_state!=State::Colliding)
        {
          //Face X
          xFaceVolume+=MathFunctions::calcCubicBSpline_Integ(0, iIndexIncrement, jIndexIncrement, kIndexIncrement);

          //Face Y
          yFaceVolume+=MathFunctions::calcCubicBSpline_Integ(1, iIndexIncrement, jIndexIncrement, kIndexIncrement);

          //Face Z
          zFaceVolume+=MathFunctions::calcCubicBSpline_Integ(2, iIndexIncrement, jIndexIncrement, kIndexIncrement);
        }

      }
    }
  }

  //Multiply interpolation integration sum with cell volume
  float cellVolume=pow(m_cellSize, 3);
  xFaceVolume*=cellVolume;
  yFaceVolume*=cellVolume;
  zFaceVolume*=cellVolume;

  //Set densities to faces
  m_cellFacesX[_cellIndex]->m_density=(m_cellFacesX[_cellIndex]->m_mass/xFaceVolume);
  m_cellFacesY[_cellIndex]->m_density=(m_cellFacesY[_cellIndex]->m_mass/yFaceVolume);
  m_cellFacesZ[_cellIndex]->m_density=(m_cellFacesZ[_cellIndex]->m_mass/zFaceVolume);

}

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcFaceVolume(int _iIndex, int _jIndex, int _kIndex, int _cellFaceDirection)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Set up to loop over neighbour cells in +-2 in each direction
    Verify that this does not include non-excisting cells

  Loop over neighbour cells
    For each cell, calculate cubicBSpline_Integ
    Add to sum cubicBSpline_Integ

  Clamp to 1, since maximum volume should be cellSize^3??

  Multiply with cellSize^3

  Return result

  ----------------------------------------------------------------------------------------------------------------
  */

  //Set storage for volume
  float faceVolume=0.0;

  //Set min and max increments for ijk
  int iMinIncrement=-2;
  int jMinIncrement=-2;
  int kMinIncrement=-2;
  int iMaxIncrement=2;
  int jMaxIncrement=2;
  int kMaxIncrement=2;

  //Make sure doesn't loop over outside cells.
  //Sort of makes the non-existing cells seem like collision cells as will add nothing to volume
  // i direction
  if (_iIndex==0)
  {
    iMinIncrement=0;
  }
  if (_iIndex==1)
  {
    iMinIncrement=-1;
  }
  if (_iIndex==(m_noCells-2))
  {
    iMaxIncrement=1;
  }
  if (_iIndex==(m_noCells-1))
  {
    iMaxIncrement=0;
  }

  //j direction
  if (_jIndex==0)
  {
    jMinIncrement=0;
  }
  if (_jIndex==1)
  {
    jMinIncrement=-1;
  }
  if (_jIndex==(m_noCells-2))
  {
    jMaxIncrement=1;
  }
  if (_jIndex==(m_noCells-1))
  {
    jMaxIncrement=0;
  }

  // k direction
  if (_kIndex==0)
  {
    kMinIncrement=0;
  }
  if (_kIndex==1)
  {
    kMinIncrement=-1;
  }
  if (_kIndex==(m_noCells-2))
  {
    kMaxIncrement=1;
  }
  if (_kIndex==(m_noCells-1))
  {
    kMaxIncrement=0;
  }


  //Loop over neighbouring cells that will give cubicBSpline_Integ!=0
  for (int kIndexIncrement=kMinIncrement; kIndexIncrement<(kMaxIncrement+1); kIndexIncrement++)
  {
    for (int jIndexIncrement=jMinIncrement; jIndexIncrement<(jMaxIncrement+1); jIndexIncrement++)
    {
      for (int iIndexIncrement=iMinIncrement; iIndexIncrement<(iMaxIncrement+1); iIndexIncrement++)
      {

        //Get index of incremented cell
        int cellIndexIncremented=MathFunctions::getVectorIndex((_iIndex+iIndexIncrement), (_jIndex+jIndexIncrement), (_kIndex+kIndexIncrement), m_noCells);

        //Check if cell is colliding. If colliding then will add nothing to cell face volume
        if (m_cellCentres[cellIndexIncremented]->m_state!=State::Colliding)
        {
          //Check which face direction is calculated for
          //Face X
          if (_cellFaceDirection==0)
          {
            faceVolume+=MathFunctions::calcCubicBSpline_Integ(0, iIndexIncrement, jIndexIncrement, kIndexIncrement);
          }
          //Face Y
          if (_cellFaceDirection==1)
          {
            faceVolume+=MathFunctions::calcCubicBSpline_Integ(1, iIndexIncrement, jIndexIncrement, kIndexIncrement);
          }
          //Face Z
          if (_cellFaceDirection==2)
          {
            faceVolume+=MathFunctions::calcCubicBSpline_Integ(2, iIndexIncrement, jIndexIncrement, kIndexIncrement);
          }
        }

      }
    }
  }

  //Clamp volume to 1
  ///Not sure if this should be done
//  faceVolume=std::min<float>(faceVolume, 1.0);

  //Multiply interpolation integration sum with cell volume
  float cellVolume=pow(m_cellSize, 3);
  faceVolume*=cellVolume;

  //Return face volume
  return faceVolume;



}
