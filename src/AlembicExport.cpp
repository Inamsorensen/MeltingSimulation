#include "AlembicExport.h"

namespace AbcG=Alembic::AbcGeom;


AlembicExport::AlembicExport(std::string _fileName)
{
  //Create alembic geometry output archive called explosion.abc
//  m_archive.reset(new AbcG::OArchive(Alembic::AbcCoreOgawa::WriteArchive(),"explosion.abc"));
  m_archive.reset(new AbcG::OArchive(Alembic::AbcCoreHDF5::WriteArchive(), _fileName));

  //Set up timesampling to be for 25fps
  AbcG::TimeSampling timeSamp(1.0f/25.0f,0.0f);

  //Get top of archive
  AbcG::OObject topObj(*m_archive.get(), AbcG::kTop);
//  AbcG::OObject archiveTop=*m_archive->getTop();

  //Add in time sampling
  Alembic::Util::uint32_t timeSampAdd=topObj.getArchive().addTimeSampling(timeSamp);


  //Reset particle position output
  m_pointSamples.reset(new AbcG::OPoints(topObj,"particles", timeSampAdd));
//  m_pointSamples.reset(new AbcG::OPoints(topObj, "particles"));

  ///To do: Add temperature property
//  m_temperatures.reset(new AbcG::OFloatArrayProperty(m_pointSamples->getSchema(), "temperature", timeSampAdd));

//  m_temperatures.reset(new AbcG::OFloatArrayProperty(m_pointSamples->getSchema().getUserProperties(), "temperature", timeSampAdd));
//  AbcG::OCompoundProperty compProp=m_pointSamples->getProperties();
//  AbcG::OCompoundProperty compProp2=m_pointSamples->getSchema().getUserProperties();
//  int noProp=compProp.getNumProperties();
//  m_temperatures.reset(new AbcG::OFloatArrayProperty(compProp, "temperature", timeSampAdd));
//  m_temperatures.reset(new AbcG::OFloatArrayProperty(compProp2, "temperature", timeSampAdd));
//  int noProp2=compProp.getNumProperties();
//  std::cout<<"test\n";
//  AbcG::OGeomBaseSchema_Points pointSchema=AbcG::OGeomBaseSchema_Points();;
//  AbcG::OCompoundProperty tempProp=pointSchema.getUserProperties();

}

void AlembicExport::exportFrame(std::vector<Imath::V3f> &_positions, std::vector<Alembic::Util::uint64_t> &_id)
{
  ///To do: Set temperature properties as well

  //Set up Imath vectors for position and id data
//  std::vector<Imath::V3f> positions;
//  std::vector<Alembic::Util::uint64_t> id;
//  std::vector<Alembic::Util::float32_t> temperature;

//  unsigned int noParticles=_noParticles;

//  for (unsigned int i=0; i<noParticles; i++)
//  {
//    float posX=_particlePositions->at(i).m_x;
//    float posY=_particlePositions->at(i).m_y;
//    float posZ=_particlePositions->at(i).m_z;
//    positions.push_back(Imath::V3f(posX, posY, posZ));
//    id.push_back(i);
//    temperature.push_back(_particleTemperatures->at(i));

//  }

  //Need to create samples from data
  AbcG::V3fArraySample samplePos(_positions);
  AbcG::UInt64ArraySample sampleId(_id);
//  AbcG::FloatArraySample sampleTemp(temperature);
  AbcG::OPointsSchema::Sample pointsSample(samplePos, sampleId);


  m_pointSamples->getSchema().set(pointsSample);
//  m_temperatures->set(AbcG::FloatArraySample(temperature));
//  m_temperatures->set(temperature);
//  m_temperatures->set(sampleTemp);

}
