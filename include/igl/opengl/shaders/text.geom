std::string text_geom_shader = R"(
    
    #version 150 core

    // http://ogldev.atspace.co.uk/www/tutorial27/tutorial27.html

    layout(points) in;
    layout(triangle_strip, max_vertices = 4) out;

    out vec2 gTexCoord;

    uniform mat4 view;
    uniform mat4 proj;

    uniform vec2 CellSize;
    uniform vec2 CellOffset;
    uniform vec2 RenderSize;
    uniform vec2 RenderOrigin;
    uniform float TextShiftFactor;

    in int vPosition[1];
    in int vCharacter[1];
    in float vOffset[1];

    void main()
    {
        // Code taken from https://prideout.net/strings-inside-vertex-buffers

        // Determine the final quad's position and size:
        vec4 P = gl_in[0].gl_Position + vec4( vOffset[0]*TextShiftFactor, 0.0, 0.0, 0.0 ); // 0.04
        vec4 U = vec4(1, 0, 0, 0) * RenderSize.x; // 1.0
        vec4 V = vec4(0, 1, 0, 0) * RenderSize.y; // 1.0

        // Determine the texture coordinates:
        int letter = vCharacter[0]; // used to be the character
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