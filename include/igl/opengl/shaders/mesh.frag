std::string text_frag_shader = R"(
    
    #version 150 

    // ---------- Rendering Basic Mesh ----------
    uniform vec4 fixed_color;
    uniform mat4 view;
    uniform mat4 proj;
    in VS_OUT {
        vec3 position_eye;
        vec3 normal_eye;
        vec4 Ksi;
        vec4 Kdi;
        vec4 Kai;
        vec2 texcoordi;
    } fs_out;
    uniform vec3 light_position_eye;
    vec3 Ls = vec3 (1, 1, 1);
    vec3 Ld = vec3 (1, 1, 1);
    vec3 La = vec3 (1, 1, 1);
    uniform sampler2D tex;
    uniform float specular_exponent;
    uniform float lighting_factor;
    uniform float texture_factor;

    // ---------- Out ----------
    out vec4 outColor;

    void main()
    {

        // ---------- Text Rendering ----------


        // ---------- Rendering Basic Mesh ----------

        vec3 Ia = 0.2 *La * vec3(fs_out.Kai);    // ambient intensity

        vec3 vector_to_light_eye = light_position_eye - fs_out.position_eye;
        vec3 direction_to_light_eye = normalize (vector_to_light_eye);
        float dot_prod = dot (direction_to_light_eye, normalize(fs_out.normal_eye));
        float clamped_dot_prod = max (dot_prod, 0.0);
        vec3 Id = Ld * vec3(fs_out.Kdi) * clamped_dot_prod;    // Diffuse intensity

        vec3 reflection_eye = reflect (-direction_to_light_eye, normalize(fs_out.normal_eye));
        vec3 surface_to_viewer_eye = normalize (-fs_out.position_eye);
        float dot_prod_specular = dot (reflection_eye, surface_to_viewer_eye);
        dot_prod_specular = float(abs(dot_prod)==dot_prod) * max (dot_prod_specular, 0.0);
        float specular_factor = pow (dot_prod_specular, specular_exponent);
        vec3 Is = Ls * vec3(fs_out.Ksi) * specular_factor;    // specular intensity

        vec4 color = vec4(lighting_factor * (Is + Id) + Ia + (1.0-lighting_factor) * 
                    vec3(fs_out.Kdi),(fs_out.Kai.a+fs_out.Ksi.a+fs_out.Kdi.a)/3);
        outColor = mix(vec4(1,1,1,1), texture(tex, fs_out.texcoordi), texture_factor) * color;
        if (fixed_color != vec4(0.0)) outColor = fixed_color;
    }

)";