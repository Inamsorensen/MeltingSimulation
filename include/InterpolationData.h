#ifndef INTERPOLATIONDATA
#define INTERPOLATIONDATA

#include <eigen3/Eigen/Core>

#include "Particle.h"

struct InterpolationData
{
  Particle* m_particle;
  float m_cubicBSpline;
  Eigen::Vector3f m_cubicBSpline_Diff;
  float m_cubicBSpline_Integ;
  float m_tightQuadStencil;
  Eigen::Vector3f m_tightQuadStencil_Diff;
};

#endif // INTERPOLATIONDATA

