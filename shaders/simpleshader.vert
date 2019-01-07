#version 330 core


layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec3 color;
layout (location = 3) in vec2 UVcoord;

out vec3 fragNormal;
out vec3 fragColor;
out vec2 UV;

uniform mat4 projection, modelview;
uniform mat3 normalMatrix;


void main()
{
    UV = UVcoord;
    fragColor = color;
    gl_Position = projection * modelview * vec4(position, 1.0);
    fragNormal = normalMatrix * normal;
}

