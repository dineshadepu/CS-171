varying vec3 normal, vec;
uniform int nLights;

void main()
{
    // transform normals to camera space and normalize
    normal = normalize(gl_NormalMatrix * gl_Normal);

    // using the transformation matrix, transform points to camera space 
    vec = vec3(gl_ModelViewMatrix * gl_Vertex); 

    gl_Position = ftransform();
}
