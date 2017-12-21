varying vec3 normal, vec;
uniform int nLights;

void main()
{ 
    // Do Phong shading based on algorithm from hw 2

    vec4 diffuse_sum = vec4(0., 0., 0., 0.);
    vec4 specular_sum = vec4(0., 0., 0., 0.);
    vec3 e_direction = normalize(-vec);

    for (int i = 0; i < nLights; i++)
    {
        vec3 l_direction = normalize(vec3(gl_LightSource[i].position)-vec);

        vec4 diffuse = gl_LightSource[i].diffuse * 
        max(dot(normal,l_direction), 0.0);
        
        vec4 specular = gl_LightSource[i].specular * 
        pow(max(dot(normal, normalize(e_direction + l_direction)), 0.0), 
            gl_FrontMaterial.shininess);

        diffuse_sum += diffuse;
        specular_sum += specular;
    }
    vec4 ambient = gl_FrontMaterial.ambient;

    vec4 color = clamp(ambient + gl_FrontMaterial.diffuse*diffuse_sum
        + gl_FrontMaterial.specular*specular_sum, 0.0, 1.0);
    gl_FragColor = color;
}