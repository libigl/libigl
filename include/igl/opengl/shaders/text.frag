std::string text_frag_shader = R"(

    #version 330

    out vec4 outColor;
    in vec2 gTexCoord;

    uniform sampler2D Sampler;

    void main()
    {
        float A = texture(Sampler, gTexCoord).r;
        outColor = vec4(1.0, 0.0, 0.0, A);
    }
)";