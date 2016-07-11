#include "ReadGeo.h"

//----------------------------------------------------------------------------------------------------------------------

ReadGeo::ReadGeo(std::string _fileName)
{
  /// @brief Opens the file for reading

  //Open file
  m_file.open(_fileName);

  if (!m_file.is_open())
  {
    std::cout<<"Failed to open file "<<_fileName<<"\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::cout<<"Opening file for reading.\n";
  }


}

//----------------------------------------------------------------------------------------------------------------------

ReadGeo::~ReadGeo()
{
  /// @brief Closes file

  if (m_file.is_open())
  {
    std::cout<<"Closing file for reading.\n";
    m_file.close();
  }
}

//----------------------------------------------------------------------------------------------------------------------

void ReadGeo::getPointPositions(int o_noPoints, std::vector<Eigen::Vector3f> &o_positionData)
{
  /// @brief Reads in number of points and position data. Done by searching for certain words in the file and reading the values after it
  /// Also checks that it has found the same number of position data as there are points according to the file

  if (m_file.is_open())
  {
    //Words to search for
    std::string pointCount="\"pointcount\"";
    std::string attributeType="\"pointattributes\"";
    std::string dataType="\"tuples\"";
    std::string paramName="P";

    //Data that is read from file
    o_noPoints=0;

    //Return pointer to beginning of file
    m_file.clear();
    m_file.seekg(0, std::ios::beg);

    //Line to read in
    std::string line;

    //Find number of points
    while (m_file>>line)
    {
      //Find number of points
      if (line.find(pointCount)!=std::string::npos)
      {
//        std::cout<<line<<"\n";

        int start=line.find_first_of(",")+1;
        int end=line.find_last_of(",");
        std::string noPointsString;
        for (int i=start; i<end; i++)
        {
          noPointsString+=line[i];
        }
        o_noPoints=std::stoi(noPointsString);
//        std::cout<<"Number of points: "<<noPoints<<"\n";

        break;
      }
    }

    //If no points found, then return.
    if (o_noPoints==0)
    {
      std::cout<<"No points found\n";
      return;
    }

    //Find line with position data
    std::string dataString;
    getDataLine(attributeType, paramName, dataType, &dataString);

//    std::cout<<dataString<<"\n";

    //Go through remaining part of line
    std::string xPos;
    std::string yPos;
    std::string zPos;

    int readDimension=0;
    int dataStringSize=dataString.size();

    for (int i=0; i<dataStringSize; i++)
    {
      //Read one letter at a time
      char letter=dataString[i];

      //If "[" then beginning of one set of position data
      if (letter=='[')
      {
        //Set to read into x
        readDimension=0;
      }
      //If "," then between x and y or y and z
      else if (letter==',')
      {
        readDimension+=1;
      }
      //If "]" then at the end of a set of position data.
      else if (letter==']')
      {
        //Turn string to floats
        float x=std::stof(xPos);
        float y=std::stof(yPos);
        float z=std::stof(zPos);

//          std::cout<<"x:"<<x<<" y:"<<y<<" z:"<<z<<"\n";

        //Store position data
        Eigen::Vector3f position;
        position(0)=x;
        position(1)=y;
        position(2)=z;
        o_positionData.push_back(position);

        //Clear strings for new read in
        xPos.clear();
        yPos.clear();
        zPos.clear();

      }
      //If not [ or ] or , then read into position data
      else
      {
        if(readDimension==0)
        {
          xPos+=letter;
        }
        else if (readDimension==1)
        {
          yPos+=letter;
        }
        else if (readDimension==2)
        {
          zPos+=letter;
        }
      }
    }

  //Check that the data stored in pointPositions is the same as the number of points
    int positionDataSize=o_positionData.size();
    std::cout<<"Number of points: "<<o_noPoints<<"\n";
    std::cout<<"Size of position data: "<<positionDataSize<<"\n";
    if (positionDataSize==o_noPoints)
    {
      std::cout<<"Same number of points as position data\n";
    }
    else
    {
      std::cout<<"Mismatch between number of points and number of position data\n";
    }


  }
}

//----------------------------------------------------------------------------------------------------------------------

void ReadGeo::getPointParameter_Float(std::string _paramName, std::vector<float> &o_data)
{
  /// @brief Read float parameters from file. Return file read to beginning as don't know what order these will be called in.

  if (m_file.is_open())
  {
    //Set words to find
    std::string pointCount="\"pointcount\"";
    std::string attributeType="\"pointattributes\"";
    std::string dataType="\"arrays\"";
    std::string paramName=_paramName;

    //Setup for read
    std::string line;
    float noPoints=0;

    //Return pointer to beginning of file
    m_file.clear();
    m_file.seekg(0, std::ios::beg);

    //Find number of points
    while (m_file>>line)
    {
      //Find number of points
      if (line.find(pointCount)!=std::string::npos)
      {
//          std::cout<<line<<"\n";

        int start=line.find_first_of(",")+1;
        int end=line.find_last_of(",");
        std::string noPointsString;
        for (int i=start; i<end; i++)
        {
          noPointsString+=line[i];
        }
        noPoints=std::stoi(noPointsString);
//          std::cout<<"Number of points: "<<noPoints<<"\n";

        break;
      }
    }

    if (noPoints==0)
    {
      std::cout<<"No points found\n";
      return;
    }

    //Get parameter data string
    std::string dataString;
    getDataLine(attributeType, paramName, dataType, &dataString);

    //Set up for reading values
    int dataStringSize=dataString.size();
    int dataNo=0;
    std::string singleDataString;

    for (int i=0; i<dataStringSize; i++)
    {
      //Read one letter at a time
      char letter=dataString[i];

//      std::cout<<letter<<"\n";

      //If "[" then beginning of one set of position data
      if (letter=='[')
      {
        //Do nothing
      }
      //If "," then next value
      else if (letter==',')
      {
        //Set current read to parameter data
        float singleDataValue=std::stof(singleDataString);
        o_data.push_back(singleDataValue);

        //Start reading into next value
        dataNo+=1;
        singleDataString.clear();

      }
      //If "]" then at the end of parameter data, so save final entry
      else if (letter==']')
      {
        //Set current read to parameter data
        float singleDataValue=std::stof(singleDataString);
        o_data.push_back(singleDataValue);

      }
      //If not [ or ] or , then read into position data
      else
      {
        singleDataString+=letter;
      }
    }

    //Check that the data stored in o_data is the same size as the number of points
      int dataSize=o_data.size();
      std::cout<<"Number of points: "<<noPoints<<"\n";
      std::cout<<"Size of data: "<<dataSize<<"\n";
      if (dataSize==noPoints)
      {
        std::cout<<"Same number of points as "<<paramName<<" data\n";
      }
      else
      {
        std::cout<<"Mismatch between number of points and number of "<<paramName<<" data\n";
      }

  }
}

//----------------------------------------------------------------------------------------------------------------------

float ReadGeo::getSimulationParameter_Float(std::string _paramName)
{
  /// @brief Searches for "globalattributes" then paramName then "arrays" to find line with value.
  /// Single value found is then returned. In case file isn't open, then returns 10000.

  //Set return value in case file isn't open
  float value=100000;

  if (m_file.is_open())
  {
    //Set words to look for
    std::string attributeType="globalattributes";
    std::string dataType="arrays";
    std::string paramName=_paramName;

    std::string dataString;

    getDataLine(attributeType, paramName, dataType, &dataString);


    int dataStringSize=dataString.size();
    std::string valueString;

    //Add every number to the string
    if (dataStringSize!=0)
    {
      for (int i=0; i<dataStringSize; i++)
      {
        char letter=dataString[i];

        if (letter!='[' && letter!=']')
        {
          valueString+=letter;
        }
      }

      //Turn string into float
      value=std::stof(valueString);
    }
  }

  return value;
}

//----------------------------------------------------------------------------------------------------------------------

Eigen::Vector3f ReadGeo::getSimulationParameter_Vec3(std::string _paramName)
{
  /// @brief Searches for "globalattributes" then paramName then "tuples" to find line with value.
  /// Single value found is then returned. In case file isn't open, then returns [0,0,0].

  Eigen::Vector3f result;

  if (m_file.is_open())
  {

    //Set words to look for
    std::string attributeType="globalattributes";
    std::string dataType="tuples";
    std::string paramName=_paramName;

    std::string dataString;
    getDataLine(attributeType, paramName, dataType, &dataString);

//    std::cout<<dataString<<"\n";

      //Read values from line
      //Set up for reading values
      int dataStringSize=dataString.size();
      std::string x_str;
      std::string y_str;
      std::string z_str;
      int dimension=0;

      for (int i=0; i<dataStringSize; i++)
      {
        //Read one letter at a time
        char letter=dataString[i];

  //      std::cout<<letter<<"\n";

        //If "[" then beginning of one set of position data
        if (letter=='[')
        {
          //Do nothing
        }
        //If "," then next value
        else if (letter==',')
        {
          //Set dimension to one higher
          dimension+=1;

        }
        //If "]" then at the end of parameter data, so save final entry
        else if (letter==']')
        {
          //Store data
          float x=std::stof(x_str);
          float y=std::stof(y_str);
          float z=std::stof(z_str);

          result(0)=x;
          result(1)=y;
          result(2)=z;

          //Set dimension to start
          dimension=0;

        }
        //If not [ or ] or , then read into position data
        else
        {
          if (dimension==0)
          {
            x_str+=letter;
          }
          else if (dimension==1)
          {
            y_str+=letter;
          }
          else if (dimension==2)
          {
            z_str+=letter;
          }
        }
      }
  }

  return result;
}

//----------------------------------------------------------------------------------------------------------------------

void ReadGeo::getDataLine(std::string _attributeType, std::string _paramName, std::string _dataType, std::string *o_data)
{
  /// @brief Uses 3 nested while loops to find the line of data. If parameter not found, returns empty o_data

  std::string line;
  std::string paramName="\"name\",\"";
  paramName.append(_paramName);
  paramName.append("\"");

  //Set pointer to start of file
  m_file.clear();
  m_file.seekg(0, std::ios::beg);


  while (m_file>>line)
  {
    if(line.find(_attributeType)!=std::string::npos)
    {
      while (m_file>>line)
      {
        if (line.find(paramName)!=std::string::npos)
        {
          while (m_file>>line)
          {
            if (line.find(_dataType)!=std::string::npos)
            {
              break;
            }
          }
          break;
        }
      }
      break;
    }
  }

  //Check if parameter was found
  if (m_file.eof())
  {
    std::cout<<"Parameter "<<_paramName<<" was not found.\n";
  }
  else
  {
//    std::cout<<line<<"\n";

    //Erase part of the line that doesn't contain the position data
    int endErase=line.find_first_of(",")+2;
    line.erase(0,endErase);

    //Make sure that all data is on one line, if not add the next line onto it
    std::string dataString;

    while (line!="]")
    {
      dataString+=line;
      m_file>>line;
    }

    //Return data line
    o_data->append(dataString);

  }
}
