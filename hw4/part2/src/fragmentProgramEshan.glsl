#version 130

uniform sampler2D colorTex;
uniform sampler2D normalTex;
varying vec3 vec, e_direction;
varying vec3 N, T, B;
varying vec2 tex_coord;
uniform int nLights;

void main()
{ 
    // using tangent, bitangent, and normal, create matrix
    // to transform light and eye vectors to surface coords
    mat3 TBN = transpose(mat3(T, B, N));

    vec3 normal = texture2D(normalTex, tex_coord).rgb;
    // map components to [-1,1]
    normal = 2.0*normal - 1.0;
    normal = normalize(normal);

    vec4 diffuse_sum = vec4(0., 0., 0., 0.);
    vec4 specular_sum = vec4(0., 0., 0., 0.);

    // multiplier for diffuse since glFrontMaterial.diffuse
    // doesnt work well
    float mult = 1.0;

    // do diffuse and specular sum calculation
    for (int i = 0; i < nLights; i++)
    {
        // transform light to surface coords
        vec3 l_direction = normalize(TBN*(vec3(gl_LightSource[i].position)-vec));

        vec4 diffuse = gl_LightSource[i].diffuse * max(dot(normal,l_direction), 0.0);
        
        // using eye direction transformed by TBN to surface coords
        vec4 specular = gl_LightSource[i].specular * pow(max(dot(normal, normalize(TBN*e_direction + l_direction)), 0.0), gl_FrontMaterial.shininess);

        diffuse_sum += diffuse;
        specular_sum += specular;
    }
    vec4 ambient = gl_FrontMaterial.ambient;

    vec4 color = clamp(ambient + mult*diffuse_sum + gl_FrontMaterial.specular*specular_sum, 0.0, 1.0);
    gl_FragColor = texture2D(colorTex, tex_coord) * color;
}