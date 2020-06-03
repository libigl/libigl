std::string text_vert_shader = R"(
    
    #version 150

    // ---------- Text Rendering ----------
    in int words;

    // ---------- Rendering Basic Mesh ----------
    uniform mat4 view;
    uniform mat4 proj;
    uniform mat4 normal_matrix;
    in vec3 position;
    in vec3 normal;
    in vec2 texcoord;
    in vec4 Ka;
    in vec4 Kd;
    in vec4 Ks;
    out VS_OUT {
        vec3 position_eye;
        vec3 normal_eye;
        vec4 Ksi;
        vec4 Kdi;
        vec4 Kai;
        vec2 texcoordi;
    } vs_out;

    void main()
    {

        // ---------- Text Rendering ----------




        // ---------- Rendering Basic Mesh ----------

        vs_out.position_eye = vec3 (view * vec4 (position, 1.0));
        vs_out.normal_eye = vec3 (normal_matrix * vec4 (normal, 0.0));
        vs_out.normal_eye = normalize(vs_out.normal_eye);
        gl_Position = proj * vec4 (vs_out.position_eye, 1.0);

        vs_out.Ksi = Ka;
        vs_out.Kdi = Kd;
        vs_out.Kai = Ks;
        vs_out.texcoordi = texcoord;
    }

)";