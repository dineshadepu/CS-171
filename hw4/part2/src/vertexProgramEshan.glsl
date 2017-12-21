varying vec3 vec, e_direction;
varying vec3 N, T, B;
varying vec2 tex_coord;
attribute vec3 tangent;

void main()
{
    // transform and normalize point vector, normal
    // and the direction in eye space
    N = normalize(gl_NormalMatrix * gl_Normal);
    vec = vec3(gl_ModelViewMatrix * gl_Vertex); 
    e_direction = normalize(-vec);

    // calculate tangent and bitangent
    T = vec3(normalize(gl_ModelViewMatrix * vec4(tangent, 1.0)));
    B = cross(T, N);

    tex_coord = gl_MultiTexCoord0.st;
    gl_Position = ftransform();
}
