#include <QMouseEvent>
#include <QGuiApplication>

#include <ngl/Camera.h>
#include <ngl/Light.h>
#include <ngl/Material.h>
#include <ngl/NGLInit.h>
#include <ngl/VAOPrimitives.h>
#include <ngl/ShaderLib.h>
#include <ngl/Transformation.h>

#include "OpenGLWindow.h"


//----------------------------------------------------------------------------------------------------------------------
/// @brief Increment for x and y translation with mouse movement
//----------------------------------------------------------------------------------------------------------------------
const static float s_INCREMENT=0.01;
//----------------------------------------------------------------------------------------------------------------------
/// @brief Increment for wheel zoom
//----------------------------------------------------------------------------------------------------------------------
const static float s_ZOOM=0.1;

//---------------------------------------------------------------------------------------------------------------------

OpenGLWindow::OpenGLWindow()
{

  //Set translation of scene to false
  m_translateScene=false;

  //Set rotation of scene to zero
  m_rotateScene=false;
  m_rotateX=0;
  m_rotateY=0;

  //Set title of window
  setTitle("Melting");

}

//---------------------------------------------------------------------------------------------------------------------

OpenGLWindow::~OpenGLWindow()
{
  //Delete all pointers
  delete m_simulationController;
//  delete m_vao;

  std::cout<<"Shutting down NGL, removing VAOs and shaders, deleting simulation controller\n";

}

//---------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::initializeGL()
{
  /// @brief Sets up camera and shaders
  /// @todo: Set up shaders
  ///        What should the timer be set to?

  //Initialise NGL library
  ngl::NGLInit::instance();

  //Set grey background for window
  glClearColor(0.4f, 0.4f, 0.4f, 1.0f);

  //Enable depth testing for drawing
  glEnable(GL_DEPTH_TEST);

  //Enable multisampling for smoother drawing
  glEnable(GL_MULTISAMPLE);

  //Camera setup
  ngl::Vec3 from(0,0,4);
  ngl::Vec3 to(0,0,0);
  ngl::Vec3 up(0,1,0);
  m_camera.set(from,to,up);
  m_camera.setShape(60,(float)720.0/576.0,0.5,150);

  //Set shaders for bounding box and particles
  ngl::ShaderLib* shaderLib=ngl::ShaderLib::instance();

  //Phong shader for particles
  shaderLib->createShaderProgram("Phong");
  shaderLib->attachShader("PhongVertex", ngl::ShaderType::VERTEX);
  shaderLib->attachShader("PhongFragment", ngl::ShaderType::FRAGMENT);
  shaderLib->loadShaderSource("PhongVertex", "shaders/Phong.vs");
  shaderLib->loadShaderSource("PhongFragment", "shaders/Phong.fs");
  shaderLib->compileShader("PhongVertex");
  shaderLib->compileShader("PhongFragment");
  shaderLib->attachShaderToProgram("Phong", "PhongVertex");
  shaderLib->attachShaderToProgram("Phong", "PhongFragment");
  shaderLib->bindAttribute("Phong", 0, "inVert");
  shaderLib->bindAttribute("Phong", 1, "inUV");
  shaderLib->bindAttribute("Phong", 2, "inNormal");
  shaderLib->linkProgramObject("Phong");
  (*shaderLib)["Phong"]->use();
  shaderLib->setShaderParam1i("Normalize", 1);

  //Create material
  ngl::Material material(ngl::STDMAT::GOLD);
  material.loadToShader("material");

  //Create light source
  ngl::Light light(ngl::Vec3(2,2,20), ngl::Colour(1,1,1,1), ngl::Colour(1,1,1,1), ngl::LightModes::POINTLIGHT);
  ngl::Mat4 IV=m_camera.getViewMatrix();
  IV.transpose();
  light.setTransform(IV);
  light.setAttenuation(1,0,0);
  light.enable();
  light.loadToShader("light");
  ngl::VAOPrimitives *vaoPrimitives=ngl::VAOPrimitives::instance();
  vaoPrimitives->createSphere("sphere", 0.1, 10);

  //Colour shader for bounding box
  shaderLib->createShaderProgram("Colour");
  shaderLib->attachShader("ColourVertex", ngl::ShaderType::VERTEX);
  shaderLib->attachShader("ColourFragment", ngl::ShaderType::FRAGMENT);
  shaderLib->loadShaderSource("ColourVertex", "shaders/colour.vs");
  shaderLib->loadShaderSource("ColourFragment", "shaders/colour.fs");
  shaderLib->compileShader("ColourVertex");
  shaderLib->compileShader("ColourFragment");
  shaderLib->attachShaderToProgram("Colour", "ColourVertex");
  shaderLib->attachShaderToProgram("Colour", "ColourFragment");
  shaderLib->linkProgramObject("Colour");
  (*shaderLib)["Colour"]->use();

  //Set up simulation controller
  m_simulationController=SimulationController::instance();

//  m_simulationController->setRenderParameters(&m_camera, "Phong");
  m_simulationController->setRenderParameters(&m_camera, "Colour");

  //Setup VAO for bounding box
  buildVAO();

  //Need to set size of initial viewport
  glViewport(0, 0, width(), height());

  //Start timer
  startTimer(10);


}

//---------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::paintGL()
{
  /// @brief Clears window, then sets new scene transformation matrix which is used to then draw bounding box
  ///        and particles

  //Clear screen
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, m_windowWidth, m_windowHeight);

  //Get scene rotations
  ngl::Mat4 rotationX;
  ngl::Mat4 rotationY;
  rotationX.rotateX(m_rotateX);
  rotationY.rotateY(m_rotateY);
  m_transformationScene=rotationX*rotationY;
  m_transformationScene.m_30=m_scenePosition.m_x;
  m_transformationScene.m_31=m_scenePosition.m_y;
  m_transformationScene.m_32=m_scenePosition.m_z;


  //Draw bounding box
  ngl::ShaderLib* shaderLib=ngl::ShaderLib::instance();
  shaderLib->use("Colour");

  ngl::Mat4 MVP;
  ngl::Mat4 M;
  ngl::Transformation modelMatrix_BoundingBox;

  //Get grid information for bounding box
  ngl::Vec3 boundingBoxPosition;
  Eigen::Vector3f bBoxPos_Eigen=m_simulationController->getBoundingBoxPosition();
  boundingBoxPosition.m_x=bBoxPos_Eigen(0);
  boundingBoxPosition.m_y=bBoxPos_Eigen(1);
  boundingBoxPosition.m_z=bBoxPos_Eigen(2);
  float boundingBoxSize=m_simulationController->getBoundingBoxSize();

  //Make bounding box slightly bigger than grid so no "flashing" when cells are visualised as well
  float incrementGrid=0.002;
  boundingBoxPosition-=ngl::Vec3(incrementGrid, incrementGrid, incrementGrid);
  boundingBoxSize+=incrementGrid*2.0;

  //Set bounding box transformation matrix
  modelMatrix_BoundingBox.setPosition(boundingBoxPosition.m_x, boundingBoxPosition.m_y, boundingBoxPosition.m_z);
  modelMatrix_BoundingBox.setScale(boundingBoxSize, boundingBoxSize, boundingBoxSize);

  //Get and set MVP
  M=modelMatrix_BoundingBox.getMatrix()*m_transformationScene;
  MVP=M*m_camera.getVPMatrix();
  shaderLib->setShaderParamFromMat4("MVP", MVP);

  //Set colour
  ngl::Vec3 colour(0.0,0.0,0.0);
  shaderLib->setRegisteredUniformVec3("colour", colour);

  m_vao->bind();
  m_vao->draw();
  m_vao->unbind();



  //Visualise grid
//  visualiseGrid();

  //Draw particles
  m_simulationController->render(m_transformationScene);


}

//---------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::resizeGL(QResizeEvent *_event)
{
  /// @brief Gets width and height from event and resets ngl camera
  ///        Why different FOV, min and max??

  //Find window width and height from resize event
  m_windowWidth=_event->size().width()*devicePixelRatio();
  m_windowWidth=_event->size().width()*devicePixelRatio();

  //Set new camera width and height
  m_camera.setShape(45.0f, (float)width()/height(), 0.05f, 350.0f);

}

//---------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::resizeGL(int _w, int _h)
{
  /// @brief Uses input width and height to reset camera

  //Set new camera widht and height
  m_camera.setShape(45.0f, (float)_w/_h, 0.05f, 350.0f);

  //Calculate new width and height
  m_windowWidth=_w*devicePixelRatio();
  m_windowHeight=_h*devicePixelRatio();

}

//---------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::buildVAO()
{
  /// @brief: Creates Vertex Array Object of a cube using lines. Sets up indices (drawing a square at a time), vertices
  ///         and colours for each vertex. Colour currently set to red.

//  ngl::VertexArrayObject::createVOA(GL_LINES);
//  m_vao= new ngl::VertexArrayObject::createVOA(GL_LINES);
  m_vao.reset(ngl::VertexArrayObject::createVOA(GL_LINES));
  m_vao->bind();

//Line indices
  const static GLubyte indices[]=  {
                                      0,1,1,2,2,3,3,0, //top
                                      0,4,4,5,5,1,1,0, //back
                                      0,4,4,7,7,3,3,0, //left
                                      3,2,2,6,6,7,7,3, //front
                                      7,6,6,5,5,4,4,7, //bottom
                                   };
//Vertices of lines
   GLfloat vertices[] = {0,1,0,
                         1,1,0,
                         1,1,1,
                         0,1,1,
                         0,0,0,
                         1,0,0,
                         1,0,1,
                         0,0,1
                        };


//Feed data into the VAO
   m_vao->setIndexedData(24*sizeof(GLfloat),vertices[0],sizeof(indices),&indices[0],GL_UNSIGNED_BYTE,GL_STATIC_DRAW);
   m_vao->setVertexAttributePointer(0,3,GL_FLOAT,0,0);
   m_vao->setNumIndices(sizeof(indices));
   m_vao->unbind();

}

//---------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::visualiseGrid()
{
  //Draw bounding box
  ngl::ShaderLib* shaderLib=ngl::ShaderLib::instance();
  shaderLib->use("Colour");

  //Get grid information
  int noCells=m_simulationController->getNoGridCells();
  Eigen::Vector3f gridPosition=m_simulationController->getGridPosition();
  float gridCellSize=m_simulationController->getGridCellSize();

  //Loop over grid cells
  for (int kIndex=0; kIndex<noCells; kIndex++)
  {
    for (int jIndex=0; jIndex<noCells; jIndex++)
    {
      for (int iIndex=0; iIndex<noCells; iIndex++)
      {
        ngl::Mat4 MVP;
        ngl::Mat4 M;
        ngl::Transformation modelMatrix_Cell;

        //Calculate cell position
        ngl::Vec3 cellPosition;
        cellPosition.m_x=gridPosition(0) + (iIndex*gridCellSize);
        cellPosition.m_y=gridPosition(1) + (jIndex*gridCellSize);
        cellPosition.m_z=gridPosition(2) + (kIndex*gridCellSize);

        //Make cells slightly smaller than they should be for ease of visualisation
        cellPosition+=ngl::Vec3(0.01, 0.01, 0.01);
        float smallerCellSize=gridCellSize-0.02;

        //Set position and size
        modelMatrix_Cell.setPosition(cellPosition.m_x, cellPosition.m_y, cellPosition.m_z);
        modelMatrix_Cell.setScale(smallerCellSize, smallerCellSize, smallerCellSize);

        //Get and set MVP
        M=modelMatrix_Cell.getMatrix()*m_transformationScene;
        MVP=M*m_camera.getVPMatrix();
        shaderLib->setShaderParamFromMat4("MVP", MVP);

        //Get cell state for visualisation
        int cellIndex=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex, noCells);

        State cellState=m_simulationController->getGridCellState(cellIndex);

        float cellTemperature=m_simulationController->getGridCellTemperature(cellIndex);
        float ambientTemp=m_simulationController->getAmbientTemperature();
        float heatSourceTemp=m_simulationController->getHeatSourceTemperature();

        //Set colours for visualisation
        ngl::Vec3 colourRed(0.8, 0.0, 0.0);
        ngl::Vec3 colourGreen(0.0, 0.8, 0.0);
        ngl::Vec3 colourBlue(0.0, 0.0, 0.8);

        if (cellState==State::Colliding)
        {
          shaderLib->setRegisteredUniformVec3("colour", colourRed);
        }
        else if (cellState==State::Interior)
        {
          shaderLib->setRegisteredUniformVec3("colour", colourGreen);
        }
        else if (cellState==State::Empty)
        {
          shaderLib->setRegisteredUniformVec3("colour", colourBlue);
        }

//        //Visualise cell temperatures
//        if (cellTemperature==heatSourceTemp)
//        {
//          shaderLib->setRegisteredUniformVec3("colour", colourRed);
//        }
//        else if (cellTemperature==ambientTemp)
//        {
//          shaderLib->setRegisteredUniformVec3("colour", colourGreen);
//        }
//        else
//        {
//          shaderLib->setRegisteredUniformVec3("colour", colourBlue);
//        }

        //Draw
        m_vao->bind();
        m_vao->draw();
        m_vao->unbind();
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::mousePressEvent(QMouseEvent* _event)
{
  /// @brief Called when a mouse button is pressed. If left, then set to rotate, if right then translate

  if (_event->button()==Qt::LeftButton)
  {
    //Set origin of rotation to current mouse position
    m_originXRot=_event->x();
    m_originYRot=_event->y();

    //Set rotation to true
    m_rotateScene=true;
  }
  else if (_event->button()==Qt::RightButton)
  {
    //Set origin of translation to current mouse position
    m_originXTrans=_event->x();
    m_originYTrans=_event->y();

    //Set translate to true
    m_translateScene=true;
  }

}

//---------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::mouseMoveEvent(QMouseEvent* _event)
{
  /// @brief Called when mousePressEvent/mouseReleaseEvent is generated
  ///        buttons() is the button state when event was called and is different from button() which checkes which
  ///        button was pressed.

  //If rotate, then measure mouse movement and set to rotation
  if (m_rotateScene==true && _event->buttons()==Qt::LeftButton)
  {
    //Find difference between current position and origin
    int diffX=(int)(_event->x()-m_originXRot);
    int diffY=(int)(_event->y()-m_originYRot);

    //Add this to scene rotation. Opposite of mouse move direction because if move in x then rotate around y-axis and vice versa
    m_rotateX+=(float)0.5f*diffY;
    m_rotateY+=(float)0.5f*diffX;

    //Set new origin for rotation
    m_originXRot=_event->x();
    m_originYRot=_event->y();

    //Update window
    update();
  }
  //If translate, measure mouse movement and set to translation
  else if (m_translateScene==true && _event->buttons()==Qt::RightButton)
  {
    //Find difference between current position and origin
    int diffX=(int)(_event->x()-m_originXTrans);
    int diffY=(int)(_event->y()-m_originYTrans);

    //Add difference to scene position
    m_scenePosition.m_x+=s_INCREMENT*diffX;
    m_scenePosition.m_y-=s_INCREMENT*diffY;

    //Set new origin for translation
    m_originXTrans=_event->x();
    m_originYTrans=_event->y();

    //Update window
    update();
  }

}

//---------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::mouseReleaseEvent(QMouseEvent* _event)
{
  /// @brief Called when mouse button released. Sets rotate or translate to false

  if (_event->button()==Qt::LeftButton)
  {
    m_rotateScene=false;
  }
  else if (_event->button()==Qt::RightButton)
  {
    m_translateScene=false;
  }

}

//---------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::wheelEvent(QWheelEvent* _event)
{
  /// @brief Called when mouse wheel is rotated. Measures the difference in wheel position and sets this to zoom

  //Check if zoom in or out
  if (_event->delta()>0)
  {
    m_scenePosition.m_z+=s_ZOOM;
  }
  else if (_event->delta()<0)
  {
    m_scenePosition.m_z-=s_ZOOM;
  }

  //Update window
  update();

}

//---------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::keyPressEvent(QKeyEvent* _event)
{
  /// @brief Called when a keyboard key is pressed.

  switch (_event->key())
  {
  case Qt::Key_Escape: QGuiApplication::exit(EXIT_SUCCESS); break;
  case Qt::Key_F: showFullScreen(); break;
  case Qt::Key_U: m_simulationController->update(); break;
  default: break;
  }

  //Update window
  update();

}

//---------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::timerEvent(QTimerEvent* _event)
{
  /// @brief Timer should update simulation. However, need to make sure time between events and simulation time step
  ///        work well. Also future export should not happen every timer event
  /// @todo Check whether sim time step and timer work well
  ///       Export timer

  //Update simulation
  m_simulationController->update();

  //Update window
  update();

}

//---------------------------------------------------------------------------------------------------------------------
