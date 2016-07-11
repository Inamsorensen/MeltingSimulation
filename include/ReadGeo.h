#ifndef READGEO
#define READGEO

#include <iostream>
#include <fstream>
#include <vector>

#include <eigen3/Eigen/Core>

#include <ngl/Vec3.h>

//------------------------------------------------------------------------------------------------------------------------------------------------------
/// @file ReadGeo.h
/// @brief Reads data from file. Reads point positions, point parameters and overall simulation parameters.
/// @author Ina M. Sorensen
/// @version 2.0
/// @date 27.06.16
///
/// @todo
//------------------------------------------------------------------------------------------------------------------------------------------------------


class ReadGeo
{
public:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Constructor. Opens file with filename _fileName for reading.
  /// @param [in] _fileName is name of file to be read
  //----------------------------------------------------------------------------------------------------------------------
  ReadGeo(std::string _fileName);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Destructor. Closes file
  //----------------------------------------------------------------------------------------------------------------------
  ~ReadGeo();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Reads point positions from file and returns them in positionData
  /// @param [out] Pointer to vector containing position data
  //----------------------------------------------------------------------------------------------------------------------
  void getPointPositions(int o_noPoints, std::vector<Eigen::Vector3f> &o_positionData);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Reads a point parameter, ie. a different value for each point.
  /// @param [in] _paramName is the name of the parameter to be read
  /// @param [out] o_data is the vector containing the parameter values
  //----------------------------------------------------------------------------------------------------------------------
  void getPointParameter_Float(std::string _paramName, std::vector<float> &o_data);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Reads a simulation parameter, ie. one for entire file. In this case a float value
  /// @param [in] _paramName is the name of the parameter to be read
  //----------------------------------------------------------------------------------------------------------------------
  float getSimulationParameter_Float(std::string _paramName);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Reads a simulation parameter, ie. one for entire file. In this case a vec3
  /// @param [in] _paramName is the name of the parameter to be read
  //----------------------------------------------------------------------------------------------------------------------
  Eigen::Vector3f getSimulationParameter_Vec3(std::string _paramName);

private:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Variable which allows the file to be read
  //----------------------------------------------------------------------------------------------------------------------
  std::ifstream m_file;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Function to find line of data in file
  //----------------------------------------------------------------------------------------------------------------------
  void getDataLine(std::string _attributeType, std::string _paramName, std::string _dataType, std::string *o_data);

};

#endif // READGEO

