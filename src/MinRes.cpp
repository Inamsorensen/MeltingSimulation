#include "MathFunctions.h"

#include <eigen3/Eigen/Core>

#include <limits>
#include <vector>
#include <cmath>

#include <iomanip>


void MathFunctions::MinRes(const Eigen::MatrixXf &_A, const Eigen::VectorXf &_B, Eigen::VectorXf &io_x, const Eigen::MatrixXf &_preconditioner, float _shift, float _maxLoops, float _tolerance, bool _show)
{
  /* Check whether A matrix and preconditioner are symmetric and that A is not singular
  -----------------------------------------------------------------------------------------------

  -----------------------------------------------------------------------------------------------
  */

  //Check that matrices are symmetric
  Eigen::MatrixXf A_transpose=_A.transpose();

//  if (A_transpose!=_A)
//  {
//    throw std::invalid_argument("A matrix is not symmetric.");
//  }

  //Set empty matrix to compare preconditioner to
  Eigen::MatrixXf emptyMatrix;

  if (_preconditioner!=emptyMatrix)
  {
    Eigen::MatrixXf precond_transpose=_preconditioner.transpose();

    if (precond_transpose!=_preconditioner)
    {
      throw std::invalid_argument("Preconditioner is not symmetric.");
    }
  }

  //Check that A matrix is not singular, ie. that it has a determinant
  float detA=_A.determinant();
  if (detA==0)
  {
    std::cout<<"A is a singular matrix. The current MINRES might not give the correct solution.\n";
  }


  /* Set up for calculation including exit messages
  ----------------------------------------------------------------------------------------------------

  ----------------------------------------------------------------------------------------------------
  */

  //Find size of system. Ie. _B->nx1 and _A->nxn where n is systemSize
  int systemSize=_B.rows();

  float minDifference_epsilon=std::numeric_limits<float>::epsilon();

    std::vector<std::string> msg(11);
    msg[0]  = " beta1 = 0.  The exact solution is  x = 0 ";
    msg[1]  = " A solution to Ax = b was found, given tol ";
    msg[2]  = " A least-squares solution was found, given tol "; //Singular solution
    msg[3]  = " Reasonable accuracy achieved, given eps "; //Need second one for singular
    msg[4]  = " x has converged to an eigenvector ";
    msg[5]  = " acond has exceeded 0.1/eps ";
    msg[6]  = " The iteration limit was reached ";
    msg[7]  = " A  does not define a symmetric matrix ";
    msg[8]  = " M  does not define a symmetric matrix ";
    msg[9]  = " M  does not define a pos-def preconditioner ";
    msg[10] = " beta2 = 0.  If M = I, b and x are eigenvectors ";

    if(_show)
    {
      std::cout<<std::setfill('-')<<std::setw(80)<<"-"<< "\n";
      std::cout<<"|            Adapted from tminres.hpp, Stanford University, 03 Jul 2016            |\n";
      std::cout<<"|                Solution of symmetric Ax=b or (A-shift*I)x = b                 |\n";
      std::cout<<std::setfill('-')<<std::setw(80)<<"-"<<"\n";
      std::cout<<std::setfill(' ');
      std::cout<<"shift = "<< _shift << "; tolerance = " << _tolerance << "; max iterations = " << _maxLoops<<"\n";
    }

    //Set up stop variables
    int stopMessage=0;
    bool calcDone=false;
    int iterations=0;

    //Set up estimates
    //||A||
    float Anorm=0.0;
    //cond(A) is the condition number of A
    float Acond=0.0;
    //||Ar_{k}||
    float Arnorm=0.0;
    //||r_{k}||
    float rnorm=0.0;
    //||y_{k}||
    float ynorm=0.0;

    /* Set up y and v for first Lanczos vector
    ----------------------------------------------------------------------------------------
    v0=0
    r_k_2=b-(A-shift*I)x0
    y=b-(A-shift*I)x0 where x0 is the initial guess - This is unless a preconditioner is applied
    beta_1=sqrt(r_k_2*y)
    v1=y/beta_1 is set inside loop
    ----------------------------------------------------------------------------------------
    */

    //Original code uses pointers for all vectors. Not sure why they need to be pointers? Quicker?
    Eigen::VectorXf r_k_2(systemSize);
    Eigen::VectorXf y(systemSize);

    ///As long as io_x is zero this is the same as setting r2=b and beta1=norm(b)

    r_k_2=(_A)*(io_x);
    r_k_2=r_k_2-(_shift*(io_x));
    r_k_2=(_B)-r_k_2;

    if (_preconditioner!=emptyMatrix)
    {
      /// @todo work out what this does. Think has to solve system My=r1
//      M->Apply(*r1, *y);
    }
    else
    {
      y=r_k_2;
    }

    float beta_1=0.0;
    beta_1=r_k_2.dot(y);

    // Test if preconditioner is a valid matrix, ie. positive definite
    if(beta_1 < 0.0)
    {
      stopMessage = 9;
      _show = true;
      calcDone = true;
    }
    else
    {
      // If b = 0 exactly stop with x = x0 as solution found.
      if(beta_1 == 0.0)
      {
        _show = true;
        calcDone = true;
      }
      else
      {
        // Set beta_1 to ||y|| as so far has been squared of y
        beta_1 = std::sqrt(beta_1);
      }
    }


    /* Initialise quantities for calculation
    ---------------------------------------------------------------------------------------
    oldbeta is the beta from the step before, ie. beta_{k-1}

    gamma, delta and epsilon are the elements of the T_k_bar matrix
    gamma_bar and delta_bar are the gamma and delta values for k+1,
    whilst oldeps is epsilon from the previous step

    phi and phi bar are the elements of Q_{k}*beta_1*e_1

    qrnorm, tnorm2 and ynorm are used in the test for convergence

    r_k_1 is the Lanczos vector from k-1 and r_k_2 the Lanczos vector from k-2

    w, w_k_1 and w_k_2 are used to calculate x. w_k_1 being w from step k-1 and w_k_2 from step k-2

    ---------------------------------------------------------------------------------------
    */

    float oldbeta=0.0;
    float beta=beta_1;
    float alpha=0.0;

    float gamma=0.0;
    float gamma_bar=0.0;
    float gamma_max=0.0;
    float gamma_min=std::numeric_limits<float>::max();
    float delta=0.0;
    float delta_bar=0.0;
    float epsilon=0.0;
    float oldeps=0.0;

    float phi=0.0;
    float phi_bar=beta_1;

    float cosinus=-1.0;  //Same as c_k
    float sinus=0.0;     //Same as s_k

    float qrnorm=beta_1;
    float tnorm2=0.0;
    float ynorm2=0.0;

//    float rhs1=beta_1;   ///Not sure what is
//    float rhs2=0.0;      ///Not sure what is
//    float z=0.0;         ///Not sure what is

    Eigen::VectorXf w(systemSize);
    Eigen::VectorXf w_k_2(systemSize);
    Eigen::VectorXf w_k_1(systemSize);
    Eigen::VectorXf r_k_1(systemSize);
    Eigen::VectorXf v(systemSize);
    w.setZero();
    w_k_2.setZero();
    w_k_1.setZero();

    r_k_1=r_k_2;



    /* Main iteration
    ----------------------------------------------------------------------------------------

    ----------------------------------------------------------------------------------------
    */
    if (!calcDone)
    {
      iterations=1;

      for (int i=0; i<_maxLoops; i++)
      {
        /* Find Lanczos vector
        ------------------------------------------------------------------
        The values calculated above: v1=(b-(A-shift*I)*io_x)/beta_1 and beta1=||v1||.
        This is unless a preconditioner is used

        for i=1,..
          w_j'=Av_j
          alpha_j=w_j'*v_j
          w_j=w_j'-alpha*v_j - beta_j*v_(j-1)
          beta_(j+1)=||w_j||
          v_(j+1)=w_j/beta_(j+1)
        end

        ------------------------------------------------------------------
        */

        //Find v from previous step
        float normaliseV=1.0/beta;
        v=normaliseV*y;

        y=(_A)*v-(_shift*v);

        //Since v0=0, only do this for steps i>0
        if (i>0)
        {
          y=y-((beta/oldbeta)*r_k_2); //Divide by oldbeta because r_k_2 not normalised
        }

        alpha=v.dot(y);
        y=y-((alpha/beta)*r_k_1);  //Divide by beta because r_k_1 not normalised

        //Set r so stores previous v vectors
        r_k_2=r_k_1;
        r_k_1=y;

        //Preconditioner
        ///@todo Need to work out what this does
        if (_preconditioner!=emptyMatrix)
        {
//          M->Apply(*r2,*y);
        }
        else
        {
          y=r_k_1; //This seems unnecessary
        }

        //Calculate new beta
        oldbeta=beta;
        beta=r_k_1.dot(y);

        //Check if preconditioner is valid
        if (beta<0)
        {
          stopMessage=9;
          break;
        }

        beta=sqrt(beta);

//        tnorm2+=alpha*alpha + oldbeta*oldbeta + beta*beta;

        //Check if beta==0, in which case an exact solution has been found
        if (i==0)
        {
          if ((beta/beta_1)<10.0*minDifference_epsilon)
          {
            stopMessage=10;
          }
        }

        /* Linear solution section
        ------------------------------------------------------------------------------
        Uses orthogonal factorisation (or QR decomposition) of tridiagonal matrix from Lanczos algorithm
           T_k=L_k_bar*Q_k
        where L_k_bar is lower triangular and Q_k is orthogonal

        Use W_k_bar and z_k_bar instead of V_k and y_k
           W_k_bar=V_k*Q_k^T
           z_k_bar=Q_k*y_k

        Then solve
           L_k_bar*z_k_bar=beta_1*e_1
           x_k=W_k_bar*z_k_bar

        Solve yk=min||beta_1*e_1 - T_{k}_bar*y||_{2} which is a least-square problem
        Need to find Q_{k}T_{k}_bar and Q_{k}beta_1*e_1 elements

        [ delta_k       epsilon_{k+1} ] = [c_k  s_k]  [delta_bar_k     0     ]
        [gamma_bar_k   delta_bar_{k+1}]   [s_k -c_k]  [  alpha_k   beta_{k+1}]

            [phi_k]     = [c_k*phi_bar_k]
        [phi_bar_{k+1}] = [s_k*phi_bar_k]

        ------------------------------------------------------------------------------
        */

        oldeps=epsilon; //Not sure if needed
        delta=cosinus*delta_bar + sinus*alpha;
        gamma_bar=sinus*delta_bar - cosinus*alpha;
        epsilon=sinus*beta;
        delta_bar=(-cosinus*beta);

        //Compute Arnorm ||Ar_{k-1}||
        float root=sqrt(gamma_bar*gamma_bar + delta_bar*delta_bar);
        Arnorm=phi_bar*root;

        //Compute c_k and s_k of Q_k for next step
        gamma=sqrt(gamma_bar*gamma_bar + beta*beta);
        gamma=std::max(gamma, minDifference_epsilon); //In case gamma is close to zero
        cosinus=gamma_bar/gamma;
        sinus=beta/gamma;

        //Compute phi and phi_bar
        phi=cosinus*phi_bar;
        phi_bar=sinus*phi_bar;


        /* Calculate new x
        ---------------------------------------------------------------------------------
        x_{k}=x_{k-1}+phi_{k}*w_{k}

        where w_{k}=(v_{k} - delta_{k}*w_{k-1} - epsilon_k*w_{k-2})/gamma_k

        ---------------------------------------------------------------------------------
        */

        float denom=1.0/gamma;

        //Save previous w vectors
        w_k_2=w_k_1;
        w_k_1=w;

        w=denom*(v-(oldeps*w_k_2)-(delta*w_k_1));

        io_x=io_x+(phi*w);


        gamma_max=std::max(gamma_max, gamma);
        gamma_min=std::min(gamma_min, gamma);


        /* Estimate norms
        ------------------------------------------------------------------------------------

        ------------------------------------------------------------------------------------
        */

        tnorm2+=alpha*alpha + oldbeta*oldbeta + beta*beta;
        //The above was moved further down as don't think any variables changes

        Anorm=sqrt(tnorm2);
        ynorm2=io_x.dot(io_x);
        ynorm=sqrt(ynorm2);

        float epsilon_A=Anorm*minDifference_epsilon;
        float epsilon_x=epsilon_A*ynorm;

        qrnorm=phi_bar;
        rnorm=qrnorm;

        float test1=0.0;
        float test2=0.0;
        test1=rnorm/(Anorm*ynorm);  //||r||/(||A|| ||x||)
        test2=root/Anorm;           //||A r_(k-1)||/(||A|| ||r_(k-1)||)  Unsure whether this should be Arnorm instead of Anorm

        //Estimate cond(A)
        Acond=gamma_max/gamma_min;


        /* Check if stopping criteria reached
        ------------------------------------------------------------------------------------

        ------------------------------------------------------------------------------------
        */

        if (stopMessage==0)
        {
          float t1=1.0+test1;
          float t2=1.0+test2;

          if (t2<=1.0)
          {
            stopMessage=2;
          }
          if (t1<=1.0)
          {
            stopMessage=1;
          }

          if (i>=(_maxLoops-1))
          {
            stopMessage=6;
          }

          if (Acond>=(0.1/minDifference_epsilon))
          {
            stopMessage=4;
          }

          if (epsilon_x>=beta_1)
          {
            stopMessage=3;
          }

          if (test2<=_tolerance)
          {
            stopMessage=2;
          }

          if (test1 <= _tolerance)
          {
            stopMessage=1;
          }

        }

        //Increase iteration count
        iterations+=1;

        if (stopMessage!=0)
        {
          break;
        }

      }


      // Display final status
      if(_show)
      {
        std::cout << std::setfill('-') << std::setw(80) << "-" << "\n";
        std::cout << msg[stopMessage] << "\n";
        std::cout << " Number of iterations: " << iterations << "\n";
        std::cout << " Anorm = " << Anorm << "\t Acond = " << Acond << "\n";
        std::cout << " rnorm = " << rnorm << "\t ynorm = " << ynorm << "\n";
        std::cout << " Arnorm = " << Arnorm << "\n";
        std::cout << std::setfill('-') << std::setw(80) << "-" << std::endl;
        std::cout << std::setfill(' ');
      }

      calcDone=true;
    }
}
