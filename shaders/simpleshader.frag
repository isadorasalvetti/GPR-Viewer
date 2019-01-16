#version 330 core

in vec3 fragNormal;
in vec3 fragColor;
in vec2 UV;

out vec4 finalColor;

uniform bool bLighting;
uniform bool useUV;


void main()
{
    vec3 myColor = vec3(0,1,0);

    if(useUV){
     float gridDivision = 50;
     int a = int(UV.x/(1/gridDivision));
     int b = int(UV.y/(1/gridDivision));

     if (a % 2 == 0 && b % 2 != 0) myColor = vec3(UV.x, 0, UV.y);
     else if (a % 2 != 0 && b % 2 == 0) myColor = vec3(UV.x, 0, UV.y);
     else myColor = vec3(1, 1, 1);
    }
    else myColor = fragColor;

    if(bLighting) finalColor = vec4(myColor * normalize(fragNormal).z, 1);
    else finalColor = vec4(myColor, 1);
}
