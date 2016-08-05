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

  set Boundary conditions

  Solve using conjugate gradient

  Set result to cells

  ----------------------------------------------------------------------------------------------------------------
  */

  //Calculate cell face densities
  calcFaceDensities();

  //Set up matrices for linear system
  Eigen::SparseMatrix<double> A_matrix(m_totNoCells, m_totNoCells);
  Eigen::VectorXd B_vector(m_totNoCells);
  Eigen::VectorXd solution(m_totNoCells);

//  Eigen::MatrixXf A_matrix_test(totNoCells, totNoCells);
//  A_matrix_test.setZero();

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
    //Only fill in non-colliding cells
    if (m_cellCentres[cellIndex]->m_state!=State::Colliding)
    {
      //Calculate B element
      B_vector(cellIndex)=calcBComponent_projectVelocity(cellIndex);

      //Insert A elements
      calcAComponent_projectVelocity(cellIndex, A_matrix);
    }

  }

  //Solve system
  float maxLoops=30;
  float minResidual=0.00001;
  MathFunctions::conjugateGradient(A_matrix, B_vector, solution, maxLoops, minResidual);


  //Use results to calculate new velocities
#pragma omp parallel for
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    int iIndex=m_cellCentres[cellIndex]->m_iIndex;
    int jIndex=m_cellCentres[cellIndex]->m_jIndex;
    int kIndex=m_cellCentres[cellIndex]->m_kIndex;

    //Only update faces that are not colliding
    if (m_cellFacesX[cellIndex]->m_state!=State::Colliding)
    {
      //Get index for cell before current cell in x direction
      int indexCellBefore=MathFunctions::getVectorIndex(iIndex-1, jIndex, kIndex, m_noCells);

      //Calculate pressure gradient
      float pressure_ijk=solution(cellIndex);
      float pressure_i_1jk=solution(indexCellBefore);

      float pressureGradient=pressure_ijk-pressure_i_1jk;

      //Calculate constant to be multiplied with gradient
      float constant=m_dt/m_cellFacesX[cellIndex]->m_density;

      //Calculate projected velocity
      m_cellFacesX[cellIndex]->m_velocity=m_cellFacesX[cellIndex]->m_velocity - (constant*pressureGradient);
    }

    //Only update faces that are not colliding
    if (m_cellFacesY[cellIndex]->m_state!=State::Colliding)
    {
      //Get index for cell before current cell in x direction
      int indexCellBefore=MathFunctions::getVectorIndex(iIndex, jIndex-1, kIndex, m_noCells);

      //Calculate pressure gradient
      float pressure_ijk=solution(cellIndex);
      float pressure_ij_1k=solution(indexCellBefore);

      float pressureGradient=pressure_ijk-pressure_ij_1k;

      //Calculate constant to be multiplied with gradient
      float constant=m_dt/m_cellFacesY[cellIndex]->m_density;

      //Calculate projected velocity
      m_cellFacesY[cellIndex]->m_velocity=m_cellFacesY[cellIndex]->m_velocity - (constant*pressureGradient);
    }

    //Only update faces that are not colliding
    if (m_cellFacesZ[cellIndex]->m_state!=State::Colliding)
    {
      //Get index for cell before current cell in x direction
      int indexCellBefore=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex-1, m_noCells);

      //Calculate pressure gradient
      float pressure_ijk=solution(cellIndex);
      float pressure_ijk_1=solution(indexCellBefore);

      float pressureGradient=pressure_ijk-pressure_ijk_1;

      //Calculate constant to be multiplied with gradient
      float constant=m_dt/m_cellFacesZ[cellIndex]->m_density;

      //Calculate projected velocity
      m_cellFacesZ[cellIndex]->m_velocity=m_cellFacesZ[cellIndex]->m_velocity - (constant*pressureGradient);
    }

  }


}

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcBComponent_projectVelocity(int _cellIndex)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Loop over all grid cells

    For empty cells set equal to zero? Or should this be done afterwards so cells here definitely zero?

    Calculate central difference of velocity for each face and add together

    Calculate constant

    Add and insert into B vector


  ----------------------------------------------------------------------------------------------------------------
  */

    //Calculate constant
    float detDeformationGradElastic=m_cellCentres[_cellIndex]->m_detDeformationGradElastic;

    //In case cells are empty, their detDeformGrad will be zero, leading to division by zero. Hence set to one.
    if (detDeformationGradElastic==0)
    {
      detDeformationGradElastic=1.0;
    }

    float constant=(detDeformationGradElastic-1.0);
    constant*=(1.0/(m_dt*detDeformationGradElastic));
    constant*=(-1.0);

    //Get indices of current cell
    int iIndex=m_cellCentres[_cellIndex]->m_iIndex;
    int jIndex=m_cellCentres[_cellIndex]->m_jIndex;
    int kIndex=m_cellCentres[_cellIndex]->m_kIndex;

    //Get index of i+1,j,k and i,j+1,k and i,j,k+1
    int index_i1jk=MathFunctions::getVectorIndex(iIndex+1, jIndex, kIndex, m_noCells);
    int index_ij1k=MathFunctions::getVectorIndex(iIndex, jIndex+1, kIndex, m_noCells);
    int index_ijk1=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex+1, m_noCells);

    //Get face velocities for all faces surrounding cell centre
    float velocityX_forward=m_cellFacesX[index_i1jk]->m_velocity;
    float velocityX_backward=m_cellFacesX[_cellIndex]->m_velocity;
    float velocityY_forward=m_cellFacesY[index_ij1k]->m_velocity;
    float velocityY_backward=m_cellFacesY[_cellIndex]->m_velocity;
    float velocityZ_forward=m_cellFacesZ[index_ijk1]->m_velocity;
    float velocityZ_backward=m_cellFacesZ[_cellIndex]->m_velocity;

    //Calculate central gradient stencil
    float centralGradient=(1.0/m_cellSize)*((velocityX_forward-velocityX_backward)+(velocityY_forward-velocityY_backward)+(velocityZ_forward-velocityZ_backward));

    //Calculate B component
    float BComponent=constant-centralGradient;

    return BComponent;

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcAComponent_projectVelocity(int _cellIndex, Eigen::SparseMatrix<double> &o_A)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------

        Check if face colliding and cell centre interior
          Get volume for each face
        Else use volume of cell

        Add, hence must be initialised to zero, GicGic' equal to (2/px +2/py + 2/pz) for i==j and -1/px/y/z otherwise x/y/z depends which value of i and j
        Need to set to zero for certain j values since only nearest neighbours have effect

        Check if colliding as need to take 1/px/y/z away from ii and set ij to zero

     End loop

     Multiply by dt/density_{i}

    End loop

  End loop

  Calculate constant

  Add constant*I - Not sure if this is faster than another loop

  ----------------------------------------------------------------------------------------------------------------
  */

  //Get indices of current cell
  int iIndex=m_cellCentres[_cellIndex]->m_iIndex;
  int jIndex=m_cellCentres[_cellIndex]->m_jIndex;
  int kIndex=m_cellCentres[_cellIndex]->m_kIndex;

  //Get indices of surrounding cells
  int cellIndex_i1jk=MathFunctions::getVectorIndex(iIndex+1, jIndex, kIndex, m_noCells);
  int cellIndex_i_1jk=MathFunctions::getVectorIndex(iIndex-1, jIndex, kIndex, m_noCells);
  int cellIndex_ij1k=MathFunctions::getVectorIndex(iIndex, jIndex+1, kIndex, m_noCells);
  int cellIndex_ij_1k=MathFunctions::getVectorIndex(iIndex, jIndex-1, kIndex, m_noCells);
  int cellIndex_ijk1=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex+1, m_noCells);
  int cellIndex_ijk_1=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex-1, m_noCells);

  //Calculate constant=dt/cellSize^2
  float constant=(m_dt/(pow(m_cellSize,2)));

  //Get densities for cell faces
  float density_i1jk=m_cellFacesX[cellIndex_i1jk]->m_density;
  float density_i_1jk=m_cellFacesX[_cellIndex]->m_density;
  float density_ij1k=m_cellFacesY[cellIndex_ij1k]->m_density;
  float density_ij_1k=m_cellFacesY[_cellIndex]->m_density;
  float density_ijk1=m_cellFacesZ[cellIndex_ijk1]->m_density;
  float density_ijk_1=m_cellFacesZ[_cellIndex]->m_density;

  //Initialise A matrix elements to zero
  float A_ijk=0.0;
  float A_i1jk=0.0;
  float A_i_1jk=0.0;
  float A_ij1k=0.0;
  float A_ij_1k=0.0;
  float A_ijk1=0.0;
  float A_ijk_1=0.0;

  //Set up sumInvDensity
  float sumInvDensity=0.0;


  //Calculate A element for surrounding cells if they are interior, ie. not colliding or empty
  switch (m_cellCentres[cellIndex_i1jk]->m_state)
  {
  case State::Empty :
  {
    sumInvDensity+=(1.0/density_i1jk);
    break;
  }
  case State::Interior :
  {
    sumInvDensity+=(1.0/density_i1jk);
    A_i1jk=(constant/density_i1jk);
    break;
  }
  default:
  {
    break;
  }
  }

  switch (m_cellCentres[cellIndex_i_1jk]->m_state)
  {
  case State::Empty :
  {
    sumInvDensity+=(1.0/density_i_1jk);
    break;
  }
  case State::Interior :
  {
    sumInvDensity+=(1.0/density_i_1jk);
    A_i_1jk=(constant/density_i_1jk);
    break;
  }
  default:
  {
    break;
  }
  }

  switch (m_cellCentres[cellIndex_ij1k]->m_state)
  {
  case State::Empty :
  {
    sumInvDensity+=(1.0/density_ij1k);
    break;
  }
  case State::Interior :
  {
    sumInvDensity+=(1.0/density_ij1k);
    A_ij1k=(constant/density_ij1k);
    break;
  }
  default:
  {
    break;
  }
  }

  switch (m_cellCentres[cellIndex_ij_1k]->m_state)
  {
  case State::Empty :
  {
    sumInvDensity+=(1.0/density_ij_1k);
    break;
  }
  case State::Interior :
  {
    sumInvDensity+=(1.0/density_ij_1k);
    A_ij_1k=(constant/density_ij_1k);
    break;
  }
  default:
  {
    break;
  }
  }

  switch (m_cellCentres[cellIndex_ijk1]->m_state)
  {
  case State::Empty :
  {
    sumInvDensity+=(1.0/density_ijk1);
    break;
  }
  case State::Interior :
  {
    sumInvDensity+=(1.0/density_ijk1);
    A_ijk1=(constant/density_ijk1);
    break;
  }
  default:
  {
    break;
  }
  }

  switch (m_cellCentres[cellIndex_ijk_1]->m_state)
  {
  case State::Empty :
  {
    sumInvDensity+=(1.0/density_ijk_1);
    break;
  }
  case State::Interior :
  {
    sumInvDensity+=(1.0/density_ijk_1);
    A_ijk_1=(constant/density_ijk_1);
    break;
  }
  default:
  {
    break;
  }
  }

  //Only if not all the surrounding cells are zero will there be a pressure for ijk
  if (sumInvDensity!=0)
  {
    A_ijk=(-constant*sumInvDensity);
  }

  //Add pressureConstant to A_ijk
  float pressureConst;
  float detDeformGradPlastic=m_cellCentres[_cellIndex]->m_detDeformationGradPlastic;
  float detDeformGradElastic=m_cellCentres[_cellIndex]->m_detDeformationGradElastic;
  float lambdaInv=m_cellCentres[_cellIndex]->m_lameLambdaInverse;

  //In case cells are empty, their detDeformGrad will be zero, leading to division by zero. Hence set to one.
  if (detDeformGradElastic==0)
  {
    detDeformGradElastic=1.0;
  }
  if (detDeformGradPlastic==0)
  {
    detDeformGradPlastic=1.0;
  }

  pressureConst=detDeformGradPlastic/detDeformGradElastic;
  pressureConst*=lambdaInv;
  pressureConst*=(1.0/m_dt);

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

//  //Insert values into matrix
//  //TEST
//  o_A_test(_cellIndex, _cellIndex)=A_ijk;
//  o_A_test(_cellIndex, cellIndex_i1jk)=A_i1jk;
//  o_A_test(_cellIndex, cellIndex_i_1jk)=A_i_1jk;
//  o_A_test(_cellIndex, cellIndex_ij1k)=A_ij1k;
//  o_A_test(_cellIndex, cellIndex_ij_1k)=A_ij_1k;
//  o_A_test(_cellIndex, cellIndex_ijk1)=A_ijk1;
//  o_A_test(_cellIndex, cellIndex_ijk_1)=A_ijk_1;


}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcFaceDensities()
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Leave for now I think

  ----------------------------------------------------------------------------------------------------------------
  */

#pragma omp parallel for
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    //Get cell indices
    int iIndex=m_cellCentres[cellIndex]->m_iIndex;
    int jIndex=m_cellCentres[cellIndex]->m_jIndex;
    int kIndex=m_cellCentres[cellIndex]->m_kIndex;

    //Set storage for density
    float xFaceDensity=0.0;
    float yFaceDensity=0.0;
    float zFaceDensity=0.0;

    //Set min and max increments for ijk
    int iMinIncrement=-2;
    int jMinIncrement=-2;
    int kMinIncrement=-2;
    int iMaxIncrement=2;
    int jMaxIncrement=2;
    int kMaxIncrement=2;

    //Make sure doesn't loop over outside cells
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
            xFaceDensity+=MathFunctions::calcCubicBSpline_Integ(0, iIndexIncrement, jIndexIncrement, kIndexIncrement);

            //Face Y
            yFaceDensity+=MathFunctions::calcCubicBSpline_Integ(1, iIndexIncrement, jIndexIncrement, kIndexIncrement);

            //Face Z
            zFaceDensity+=MathFunctions::calcCubicBSpline_Integ(2, iIndexIncrement, jIndexIncrement, kIndexIncrement);
          }

        }
      }
    }

    //Multiply interpolation integration sum with cell volume
    float cellVolume=pow(m_cellSize, 3);
    xFaceDensity*=cellVolume;
    yFaceDensity*=cellVolume;
    zFaceDensity*=cellVolume;

    //Set densities to faces
    m_cellFacesX[cellIndex]->m_density=xFaceDensity;
    m_cellFacesY[cellIndex]->m_density=yFaceDensity;
    m_cellFacesZ[cellIndex]->m_density=zFaceDensity;

  }
}
