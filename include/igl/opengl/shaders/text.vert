std::string text_vert_shader = R"(
    
    #version 330
    
    in vec3 position;
    in float character;
    in float offset;

    uniform mat4 view;
    uniform mat4 proj;

    out int vPosition;
    out int vCharacter;
    out float vOffset;

    void main()
    {
        vCharacter = int(character);
        vOffset = offset;
        vPosition = gl_VertexID;
        gl_Position = proj * view * vec4(position, 1.0);
    }

)";
