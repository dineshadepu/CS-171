// Is lighting per-pixel or per-vertex?
uniform float mode;

// The position and normal at this vertex, from our buffers
attribute vec3 v_v;
attribute vec3 n_v;

// The vectors to be interpolated for each fragment, using barycentric
// coordinates
varying vec4 v_s;
varying vec3 n_s;

void main(void) {
    // Transform the vertex into camera space
    v_s = gl_ModelViewMatrix * vec4(v_v, 1.0);

    // Project it into NDC
    gl_Position = gl_ProjectionMatrix * v_s;

    // Transform the normal into camera space
    n_s = normalize(gl_NormalMatrix * n_v);

    // If we're in per-vertex mode, calculate the Phong lighting now
    if (mode == 0.0) {
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
        v_s = clamp(a + d + s, 0.0, 1.0);
    }
}
