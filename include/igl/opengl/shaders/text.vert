std::string text_vert_shader = R"(
    
    #version 330
    
    // ---------- Text Rendering ----------
    in vec3 position;
    in int text;

    uniform mat4 view;
    uniform mat4 proj;

    out int vPosition;
    out int vCharacter;

    void main()
    {

        // ---------- Text Rendering ----------
        vCharacter = text;
        vPosition = gl_VertexID;
        gl_Position = proj * view * vec4(position, 1.0);

    }

)";