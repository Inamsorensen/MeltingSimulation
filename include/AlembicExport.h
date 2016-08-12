#ifndef ALEMBICEXPORT
#define ALEMBICEXPORT


#include <Alembic/AbcCoreAbstract/All.h>
#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreOgawa/All.h>
#include <Alembic/AbcCoreHDF5/All.h>

#include <ngl/Vec3.h>

class AlembicExport
{
public:
  //---------------------------------------------------------------------------------
  /// @brief Alembic export ctor
  //---------------------------------------------------------------------------------
  AlembicExport(std::string _fileName);
  //---------------------------------------------------------------------------------
  /// @brief Alembic export dtor
  //---------------------------------------------------------------------------------
  ~AlembicExport(){;}
  //---------------------------------------------------------------------------------
  /// @brief Write frame to alembic
  //---------------------------------------------------------------------------------
  void exportFrame(std::vector<Imath::V3f> &_positions, std::vector<Alembic::Util::uint64_t> &_id);

private:
  //---------------------------------------------------------------------------------
  /// @brief Archive for the alembic data
  //---------------------------------------------------------------------------------
  std::unique_ptr <Alembic::AbcGeom::OArchive> m_archive;
  //---------------------------------------------------------------------------------
  /// @brief Point samples given by simulation
  //---------------------------------------------------------------------------------
  std::unique_ptr <Alembic::AbcGeom::OPoints> m_pointSamples;
  //---------------------------------------------------------------------------------
  /// @brief Array of temperatures for each point
  //---------------------------------------------------------------------------------
  std::unique_ptr <Alembic::AbcGeom::OFloatArrayProperty> m_temperatures;
};

#endif // ALEMBICEXPORT

