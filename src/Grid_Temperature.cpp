#include "Grid.h"

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcTemperature()
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------

  Set up A and B and T

  Store current temp as previous

  Set B

  Set A

  Solve conjugate gradient

  Store result
  ----------------------------------------------------------------------------------------------------------------
  */

  //Set up matrices for linear system
  Eigen::SparseMatrix<double> A_matrix(m_totNoCells, m_totNoCells);
  Eigen::VectorXd B_vector(m_totNoCells);
  Eigen::VectorXd solution(m_totNoCells);

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
    //Store current temp as previous
    m_cellCentres[cellIndex]->m_previousTemperature=m_cellCentres[cellIndex]->m_temperature;

    //Only update interior cells
    if (m_cellCentres[cellIndex]->m_state==State::Interior)
    {
      //Get ijk
      int iIndex=m_cellCentres[cellIndex]->m_iIndex;
      int jIndex=m_cellCentres[cellIndex]->m_jIndex;
      int kIndex=m_cellCentres[cellIndex]->m_kIndex;

      //Calculate B element
      B_vector(cellIndex)=calcBComponent_temperature(cellIndex);

      //Insert A elements
      calcAComponent_temperature(cellIndex, iIndex, jIndex, kIndex, A_matrix);
    }

  }

  //Solve system
  float maxLoops=3000;
  float minResidual=0.00001;
  MathFunctions::conjugateGradient(A_matrix, B_vector, solution, maxLoops, minResidual);


  //Update temperature
#pragma omp parallel for
  for (int cellIndex=0; cellIndex<m_totNoCells; cellIndex++)
  {
    //Only update interior cells
    if (m_cellCentres[cellIndex]->m_state==State::Interior)
    {
      m_cellCentres[cellIndex]->m_temperature=solution(cellIndex);

      //Need to divide by mass again
//      float mass=m_cellCentres[cellIndex]->m_mass;
//      m_cellCentres[cellIndex]->m_temperature=solution(cellIndex)/mass;
    }
  }

}

//----------------------------------------------------------------------------------------------------------------------

float Grid::calcBComponent_temperature(int _cellIndex)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Loop over grid cells
    Set B[i]=Tc^{n}

  ----------------------------------------------------------------------------------------------------------------
  */

  float bComponent=0.0;

  float mass=m_cellCentres[_cellIndex]->m_mass;
  float heatCapacity=m_cellCentres[_cellIndex]->m_heatCapacity;
  float temperature=m_cellCentres[_cellIndex]->m_temperature;

//  bComponent=temperature;
  bComponent=mass*heatCapacity*temperature;
  return bComponent;

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcAComponent_temperature(int _cellIndex, int _iIndex, int _jIndex, int _kIndex, Eigen::SparseMatrix<double> &o_A)
{
  /* Outline
  ----------------------------------------------------------------------------------------------------------------
  Calculate constant
    constant=(dt*cellSize^3)/(mass*heatCapacity)

  Set X,Y,Z values of A[ijk, ijk] to -2 initially

  Calculate surrounding A components
    constant*heatConductivity_in_direction      if cell is interior
    0 and remove one from X/Y/Z of A[ijk,ijk]   if cell is empty to enforce Neumann boundary condition
    0                                           if cell is colliding leaving its temperature fixed

  Calculate A[ijk,ijk]
    A[ijk,ijk]=1.0 + constant*(heatConductivityX*A_X + heatConductivityY*A_Y + heatConductivityZ*A_Z)

  Insert A components into matrix

  ----------------------------------------------------------------------------------------------------------------
  */

  //Get indices of surrounding cells
  int cellIndex_i1jk=MathFunctions::getVectorIndex(_iIndex+1, _jIndex, _kIndex, m_noCells);
  int cellIndex_i_1jk=MathFunctions::getVectorIndex(_iIndex-1, _jIndex, _kIndex, m_noCells);
  int cellIndex_ij1k=MathFunctions::getVectorIndex(_iIndex, _jIndex+1, _kIndex, m_noCells);
  int cellIndex_ij_1k=MathFunctions::getVectorIndex(_iIndex, _jIndex-1, _kIndex, m_noCells);
  int cellIndex_ijk1=MathFunctions::getVectorIndex(_iIndex, _jIndex, _kIndex+1, m_noCells);
  int cellIndex_ijk_1=MathFunctions::getVectorIndex(_iIndex, _jIndex, _kIndex-1, m_noCells);

  //Get variables required
  float volume=pow(m_cellSize,3);
  float mass=m_cellCentres[_cellIndex]->m_mass;
  float heatCapacity=m_cellCentres[_cellIndex]->m_heatCapacity;
  float heatConductivityX=m_cellFacesX[_cellIndex]->m_heatConductivity;
  float heatConductivityY=m_cellFacesY[_cellIndex]->m_heatConductivity;
  float heatConductivityZ=m_cellFacesZ[_cellIndex]->m_heatConductivity;

  //Calculate constant
//  float constant=(m_dt*volume)/(mass*heatCapacity);
  float constant=(m_dt*volume);

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
  case State::Empty :
  {
//    A_ijk_X+=1.0;
    A_ijk_X-=1.0;
    break;
  }
  case State::Interior :
  {
//    A_i1jk=(constant*heatConductivityX);
    A_i1jk=(-1.0*constant*heatConductivityX);
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
//    A_ijk_X+=1.0;
    A_ijk_X-=1.0;
    break;
  }
  case State::Interior :
  {
//    A_i_1jk=(constant*heatConductivityX);
    A_i_1jk=(-1.0*constant*heatConductivityX);
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
//    A_ijk_Y+=1.0;
    A_ijk_Y-=1.0;
    break;
  }
  case State::Interior :
  {
//    A_ij1k=(constant*heatConductivityY);
    A_ij1k=(-1.0*constant*heatConductivityY);
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
//    A_ijk_Y+=1.0;
    A_ijk_Y-=1.0;
    break;
  }
  case State::Interior :
  {
//    A_ij_1k=(constant*heatConductivityY);
    A_ij_1k=(-1.0*constant*heatConductivityY);
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
//    A_ijk_Z+=1.0;
    A_ijk_Z-=1.0;
    break;
  }
  case State::Interior :
  {
//    A_ijk1=(constant*heatConductivityZ);
    A_ijk1=(-1.0*constant*heatConductivityZ);
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
//    A_ijk_Z+=1.0;
    A_ijk_Z-=1.0;
    break;
  }
  case State::Interior :
  {
//    A_ijk_1=(constant*heatConductivityZ);
    A_ijk_1=(-1.0*constant*heatConductivityZ);
    break;
  }
  default:
  {
    break;
  }
  }

  //Calculate A_ijk from X,Y,Z contributions, density and constant
  A_ijk=((A_ijk_X*heatConductivityX)+(A_ijk_Y*heatConductivityY)+(A_ijk_Z*heatConductivityZ));
  A_ijk*=constant;

//  //Add one to A_ijk
//  A_ijk+=1.0;

  //Add mass*heatCapacity to A_ijk
  A_ijk+=(mass*heatCapacity);

  //Insert values into matrix
  //Might need to do this in main function for parallelisation
  o_A.insert(_cellIndex, _cellIndex)=A_ijk;
  o_A.insert(_cellIndex, cellIndex_i1jk)=A_i1jk;
  o_A.insert(_cellIndex, cellIndex_i_1jk)=A_i_1jk;
  o_A.insert(_cellIndex, cellIndex_ij1k)=A_ij1k;
  o_A.insert(_cellIndex, cellIndex_ij_1k)=A_ij_1k;
  o_A.insert(_cellIndex, cellIndex_ijk1)=A_ijk1;
  o_A.insert(_cellIndex, cellIndex_ijk_1)=A_ijk_1;

}

//----------------------------------------------------------------------------------------------------------------------
