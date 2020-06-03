std::string text_geom_shader = R"(
    
    #version 150 core

    layout(triangles) in;
    layout(triangle_strip, max_vertices = 5) out;

    // ---------- Rendering Basic Mesh ----------
    in VS_OUT {
        vec3 position_eye;
        vec3 normal_eye;
        vec4 Ksi;
        vec4 Kdi;
        vec4 Kai;
        vec2 texcoordi;
    } gs_in[];
    out VS_OUT {
        vec3 position_eye;
        vec3 normal_eye;
        vec4 Ksi;
        vec4 Kdi;
        vec4 Kai;
        vec2 texcoordi;
    } gs_out;

    void main()
    {

        // ---------- Text Rendering ----------


        // ---------- Rendering Basic Mesh ----------

        gs_out.position_eye = gs_in[0].position_eye;
        gs_out.normal_eye   = gs_in[0].normal_eye;
        gs_out.texcoordi    = gs_in[0].texcoordi;
        gs_out.Kai          = gs_in[0].Kai;
        gs_out.Kdi          = gs_in[0].Kdi;
        gs_out.Ksi          = gs_in[0].Ksi;
    
        gl_Position = gl_in[1].gl_Position;
        EmitVertex();
        gl_Position = gl_in[2].gl_Position;
        EmitVertex();
        gl_Position = gl_in[0].gl_Position;
        EmitVertex();
        EndPrimitive();
    }
    
)";