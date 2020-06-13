std::string text_frag_shader = R"(

    #version 330

    out vec4 outColor;
    in vec2 gTexCoord;

    uniform sampler2D Sampler;
    uniform vec3 TextColor;

    void main()
    {
        float A = texture(Sampler, gTexCoord).r;
        outColor = vec4(TextColor, A);
    }
)";
