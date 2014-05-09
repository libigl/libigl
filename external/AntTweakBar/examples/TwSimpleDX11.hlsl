//  ---------------------------------------------------------------------------
//
//  @file       TwSimpleDX11.hlsl
//  @brief      Shaders used by the TwSimpleDX11 example
//
//  ---------------------------------------------------------------------------

cbuffer Constants : register(b0)
{
    row_major float4x4 g_WorldViewProj;
	row_major float4x4 g_WorldNorm;
    float3 g_LightDir;
    float g_LightCoeff;
};

struct PSInput 
{ 
    float4 Pos : SV_POSITION; 
    float4 Color : COLOR0; 
};

PSInput MainVS(float4 pos : POSITION, float3 norm : NORMAL, float4 color : COLOR) 
{
    PSInput ps; 
    ps.Pos = mul(pos, g_WorldViewProj);
    float3 n = normalize(mul(float4(norm, 0), g_WorldNorm).xyz);
    ps.Color.rgb = color.rgb * ((1 - g_LightCoeff) + g_LightCoeff * dot(n, -g_LightDir));
    ps.Color.a  = color.a;
    return ps;
}

float4 MainPS(PSInput input) : SV_TARGET
{
    return input.Color; 
}
