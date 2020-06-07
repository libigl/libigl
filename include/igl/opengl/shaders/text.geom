std::string text_geom_shader = R"(
    
    #version 150 core

    layout(points) in;
    layout(triangle_strip, max_vertices = 4) out;

    out vec2 gTexCoord;
    uniform sampler2D Sampler;

    uniform mat4 view;
    uniform mat4 proj;

    uniform vec2 CellSize;
    uniform vec2 CellOffset;
    uniform vec2 RenderSize;
    uniform vec2 RenderOrigin;

    in int vPosition[1];
    in int vCharacter[1];
    in float vOffset[1];

    void main()
    {
        // Determine the final quad's position and size:
        // vec4 P = gl_in[0].gl_Position;


        // float x = RenderOrigin.x + float(vPosition[0]) * RenderSize.x * 2;
        // float y = RenderOrigin.y;
        // vec4 P = vec4(x, y, 0, 1);
        
        
        // vec4 P = gl_in[0].gl_Position;
        vec4 P = gl_in[0].gl_Position + vec4( vOffset[0]*0.04, 0.0, 0.0, 0.0 );
        // P = proj * view * P;


        vec4 U = vec4(1, 0, 0, 0) * RenderSize.x;
        vec4 V = vec4(0, 1, 0, 0) * RenderSize.y;

        // Determine the texture coordinates:
        int letter = vCharacter[0]; // used to be the character

        // int letter = 97; // used to be the character

        letter = clamp(letter - 32, 0, 96);
        int row = letter / 16 + 1;
        int col = letter % 16;
        float S0 = CellOffset.x + CellSize.x * col;
        float T0 = CellOffset.y + 1 - CellSize.y * row;
        float S1 = S0 + CellSize.x - CellOffset.x;
        float T1 = T0 + CellSize.y;

        // Output the quad's vertices:
        gTexCoord = vec2(S0, T1); gl_Position = P - U - V; EmitVertex();
        gTexCoord = vec2(S1, T1); gl_Position = P + U - V; EmitVertex();
        gTexCoord = vec2(S0, T0); gl_Position = P - U + V; EmitVertex();
        gTexCoord = vec2(S1, T0); gl_Position = P + U + V; EmitVertex();
        EndPrimitive();
    }
    
)";