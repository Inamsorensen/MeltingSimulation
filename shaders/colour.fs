#version 150
/// @file Colour.fs
/// @brief a basic unshaded solid colour shader
/// @brief the colour to shade draw with
in vec3 colourOut;
uniform vec4 Colour;
out vec4 fragColour;

void main ()
{
  fragColour = vec4(colourOut,1);
}

