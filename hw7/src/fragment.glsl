// Is lighting per-pixel or per-vertex?
uniform float mode;

// The values calculated at each vertex, interpolated for this fragment
varying vec4 v_s;
varying vec3 n_s;

void main (void) {
    // If per-pixel, apply the Phong shading algorithm to each fragment (pixel)
    // using the interpolated values
    if (mode == 1.0) {
        vec4 a = vec4(0.0, 0.0, 0.0, 0.0);
        vec4 d = vec4(0.0, 0.0, 0.0, 0.0);
        vec4 s = vec4(0.0, 0.0, 0.0, 0.0);
        vec3 l, e, r;

        // Loop through all the lights in the scene
        for (int i = 0; i < gl_MaxLights; i++) {
            // First we get our l, e, and r vectors (see HW 2 lecture notes)
            l = normalize(vec3(gl_LightSource[i].position - v_s));   
            e = normalize(-vec3(v_s));  
            r = normalize(reflect(-l, n_s));
         
            // Calculate ambient term  
            a += gl_FrontLightProduct[i].ambient;    

            // Calculate diffuse term  
            d += gl_FrontLightProduct[i].diffuse * max(dot(n_s, l), 0.0);    
           
            // Calculate specular term 
            s += gl_FrontLightProduct[i].specular *
               pow(max(dot(r, e), 0.0), gl_FrontMaterial.shininess);
        }
        gl_FragColor = clamp(a + d + s, 0.0, 1.0);
    }
    // If per-vertex, interpolated color is in v_s;
    else {
        gl_FragColor = v_s;
    }
}