#version 330 core


in vec3 fragNormal;
in vec3 fragColor;

out vec4 finalColor;


uniform bool bLighting;


void main()
{
    if(bLighting)
            finalColor = vec4(fragColor * normalize(fragNormal).z, 1);
    else
            finalColor = vec4(fragColor, 1);
}
