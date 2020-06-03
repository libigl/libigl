std::string text_vert_shader = R"(
    
    #version 330
    
    // ---------- Text Rendering ----------
    layout(location = 0) in vec3  position;
    layout(location = 1) in vec3  color;
    layout(location = 2) in float text;

    uniform mat4 view;
    uniform mat4 proj;

    out int vPosition;

    void main()
    {

        // mat4 aMat4 = mat4(0.7, 0.0, 0.0, 0.0,  // 1. column
        //                   0.0, 0.7, 0.0, 0.0,  // 2. column
        //                   0.0, 0.0, 0.7, 0.0,  // 3. column
        //                   0.0, 0.0, 0.0, 0.7); // 4. column

        // ---------- Text Rendering ----------
        vPosition = gl_VertexID;
        gl_Position =  proj * view * vec4(position, 1.0);


    }

)";