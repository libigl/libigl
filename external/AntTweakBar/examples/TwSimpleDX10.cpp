//  ---------------------------------------------------------------------------
//
//  @file       TwSimpleDX10.cpp
//  @brief      A simple example that uses AntTweakBar with DirectX10.
//              This example draws a triangle and allow the user to tweak its
//              vertex positions and colors.
//
//              AntTweakBar: http://anttweakbar.sourceforge.net/doc
//              DirectX:     http://msdn.microsoft.com/directx
//  
//  @author     Philippe Decaudin
//
//  ---------------------------------------------------------------------------

#include <AntTweakBar.h>
#include <cmath>

#include "d3d10vs2003.h" // workaround to include D3D10.h with VS2003
#include <d3d10.h>


// D3D interface pointers
ID3D10Device *          g_D3DDevice = NULL;
IDXGISwapChain *        g_SwapChain = NULL;
DXGI_SWAP_CHAIN_DESC    g_SwapChainDesc;
ID3D10RenderTargetView *g_RenderTargetView = NULL;
ID3D10Effect *          g_Effect = NULL;
ID3D10EffectTechnique * g_Technique = NULL;
ID3D10InputLayout *     g_VertexLayout = NULL;
ID3D10Buffer *          g_VertexBuffer = NULL;
ID3D10BlendState *      g_BlendState = NULL;
ID3D10RasterizerState * g_RasterState = NULL;

const int NB_VERTS = 3;

 
// Forward declarations
HRESULT             InitDevice(HWND wnd);
HRESULT             InitScene();
void                Cleanup();
LRESULT CALLBACK    MessageProc(HWND, UINT, WPARAM, LPARAM);
void                Render();

// Variables
int   g_Angle = 0;
float g_Scale = 1;
struct Point { float X, Y; };
Point g_Positions[NB_VERTS] = { {0.0f, 0.5f}, {0.5f, -0.5f}, {-0.5f, -0.5f} };
float g_Colors[NB_VERTS][4] = { {0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1} };


// Main
int WINAPI WinMain(HINSTANCE instance, HINSTANCE, LPSTR, int cmdShow)
{
    // Register our window class
    WNDCLASSEX wcex = { sizeof(WNDCLASSEX), CS_HREDRAW|CS_VREDRAW, MessageProc,
                        0L, 0L, instance, NULL, NULL, NULL, NULL, L"TwDX10", NULL };
    RegisterClassEx(&wcex);

    // Create a window
    RECT rc = { 0, 0, 640, 480 };
    AdjustWindowRect(&rc, WS_OVERLAPPEDWINDOW, FALSE);
    HWND wnd = CreateWindow(L"TwDX10", L"AntTweakBar simple example using DirectX10", 
                            WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 
                            rc.right-rc.left, rc.bottom-rc.top, NULL, NULL, instance, NULL);
    if( !wnd )
    {
        MessageBox(NULL, L"Cannot create window", L"Error", MB_OK|MB_ICONERROR);
        return 0;
    }

    // Initialize D3D10
    if( FAILED(InitDevice(wnd)) )
    {
        Cleanup();
        MessageBox(wnd, L"Cannot create D3D10 device", L"Error", MB_OK|MB_ICONERROR);
        return 0;
    }

    // Initialize the 3D scene
    if( FAILED(InitScene()) )
    {
        Cleanup();
        MessageBox(wnd, L"Scene initialization failed.", L"Error", MB_OK|MB_ICONERROR);
        return 0;
    }


    // Initialize AntTweakBar
    if( !TwInit(TW_DIRECT3D10, g_D3DDevice) )
    {
        Cleanup();
        MessageBoxA(wnd, TwGetLastError(), "AntTweakBar initialization failed", MB_OK|MB_ICONERROR);
        return 0;
    }

    // Create a tweak bar
    TwBar *bar = TwNewBar("TweakBar");
    TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar in a DirectX10 application.' "); // Message added to the help bar.

    // Add rotation and scale vars to bar
    TwAddVarRW(bar, "Rotation", TW_TYPE_INT32, &g_Angle, 
               " KeyIncr=r KeyDecr=R Help='Rotates the triangle (angle in degree).' ");
    TwAddVarRW(bar, "Scale", TW_TYPE_FLOAT, &g_Scale, 
               " Min=-2 Max=2 Step=0.01 KeyIncr=s KeyDecr=S Help='Scales the triangle (1=original size).' ");

    // Create a new TwType to edit 2D points: a struct that contains two floats
    TwStructMember pointMembers[] = { 
        { "X", TW_TYPE_FLOAT, offsetof(Point, X), " Min=-1 Max=1 Step=0.01 " },
        { "Y", TW_TYPE_FLOAT, offsetof(Point, Y), " Min=-1 Max=1 Step=0.01 " } };
    TwType pointType = TwDefineStruct("POINT", pointMembers, 2, sizeof(Point), NULL, NULL);

    // Add color and position of the 3 vertices of the triangle
    TwAddVarRW(bar, "Color0", TW_TYPE_COLOR4F, &g_Colors[0], " Alpha HLS Group='Vertex 0' Label=Color ");
    TwAddVarRW(bar, "Pos0", pointType, &g_Positions[0], " Group='Vertex 0' Label='Position' ");
    TwAddVarRW(bar, "Color1", TW_TYPE_COLOR4F, &g_Colors[1], " Alpha HLS Group='Vertex 1' Label=Color ");
    TwAddVarRW(bar, "Pos1", pointType, &g_Positions[1], " Group='Vertex 1' Label='Position' ");
    TwAddVarRW(bar, "Color2", TW_TYPE_COLOR4F, &g_Colors[2], " Alpha HLS Group='Vertex 2' Label=Color ");
    TwAddVarRW(bar, "Pos2", pointType, &g_Positions[2], " Group='Vertex 2' Label='Position' ");


    // Show and update the window
    ShowWindow(wnd, cmdShow);
    UpdateWindow(wnd);

    // Main message loop.
    // Passive loop: Content is repaint only when needed.
    MSG msg = {0};
    while( GetMessage(&msg, NULL, 0, 0) > 0 )
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    TwTerminate();
    Cleanup();

    return (int)msg.wParam;
}


// Create Direct3D device and swap chain
HRESULT InitDevice(HWND wnd)
{
    HRESULT hr = S_OK;
    RECT rc;
    GetClientRect(wnd, &rc);
    UINT width = rc.right - rc.left;
    UINT height = rc.bottom - rc.top;

    UINT createDeviceFlags = 0;
    #ifdef _DEBUG
        createDeviceFlags |= D3D10_CREATE_DEVICE_DEBUG;
    #endif
    ZeroMemory(&g_SwapChainDesc, sizeof(g_SwapChainDesc));
    g_SwapChainDesc.BufferCount = 1;
    g_SwapChainDesc.BufferDesc.Width = width;
    g_SwapChainDesc.BufferDesc.Height = height;
    g_SwapChainDesc.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    g_SwapChainDesc.BufferDesc.RefreshRate.Numerator = 60;
    g_SwapChainDesc.BufferDesc.RefreshRate.Denominator = 1;
    g_SwapChainDesc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    g_SwapChainDesc.OutputWindow = wnd;
    g_SwapChainDesc.SampleDesc.Count = 1;
    g_SwapChainDesc.SampleDesc.Quality = 0;
    g_SwapChainDesc.Windowed = TRUE;
    g_SwapChainDesc.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;
   
    // Try to create a hardware accelerated device first
    hr = D3D10CreateDeviceAndSwapChain(NULL, D3D10_DRIVER_TYPE_HARDWARE, NULL, createDeviceFlags, 
                                       D3D10_SDK_VERSION, &g_SwapChainDesc, &g_SwapChain, &g_D3DDevice);
    if( FAILED(hr) )
    {
        // If failed, try to create a reference device
        hr = D3D10CreateDeviceAndSwapChain(NULL, D3D10_DRIVER_TYPE_REFERENCE, NULL, createDeviceFlags,
                                           D3D10_SDK_VERSION, &g_SwapChainDesc, &g_SwapChain, &g_D3DDevice);
        if( SUCCEEDED(hr) )
            MessageBox(wnd, L"No DX10 hardware acceleration found.\nSwitching to REFERENCE driver (very slow).",
                       L"Warning", MB_OK|MB_ICONWARNING);
        else
            return hr;
    }

    // Create a render target view
    ID3D10Texture2D *buffer;
    hr = g_SwapChain->GetBuffer( 0, __uuidof( ID3D10Texture2D ), (LPVOID*)&buffer );
    if( FAILED(hr) )
        return hr;

    hr = g_D3DDevice->CreateRenderTargetView( buffer, NULL, &g_RenderTargetView );
    buffer->Release();
    if( FAILED(hr) )
        return hr;

    g_D3DDevice->OMSetRenderTargets( 1, &g_RenderTargetView, NULL );

    // Setup the viewport
    D3D10_VIEWPORT vp;
    vp.Width = width;
    vp.Height = height;
    vp.MinDepth = 0.0f;
    vp.MaxDepth = 1.0f;
    vp.TopLeftX = 0;
    vp.TopLeftY = 0;
    g_D3DDevice->RSSetViewports( 1, &vp );

    return S_OK;
}


// Initialize the 3D objects & effects
HRESULT InitScene()
{
    HRESULT hr = S_OK;

    // Create effect
    DWORD shaderFlags = D3D10_SHADER_ENABLE_STRICTNESS;
    #if defined( DEBUG ) || defined( _DEBUG )
    //  shaderFlags |= D3D10_SHADER_DEBUG;
    #endif

    // DX10 requires shaders (no fixed function anymore). We will read the shaders from the following string.
    char shaderFX[] = " struct PSInput { float4 Pos : SV_POSITION; float4 Col : COLOR0; }; \n"
                      " PSInput VertShad(float4 pos : POSITION, float4 col : COLOR) { \n"
                      "   PSInput ps; ps.Pos=pos; ps.Col=col; return ps; } \n"
                      " float4 PixShad(PSInput input) : SV_Target { return input.Col; } \n"
                      " technique10 Render { pass P0 { \n"
                      "   SetVertexShader( CompileShader( vs_4_0, VertShad() ) ); \n"
                      "   SetGeometryShader( NULL ); \n"
                      "   SetPixelShader( CompileShader( ps_4_0, PixShad() ) ); \n"
                      " } }\n";

    ID3D10Blob *compiledFX = NULL;
    ID3D10Blob *errors = NULL;
    hr = D3D10CompileEffectFromMemory(shaderFX, strlen(shaderFX), "TestDX10", 
                                      NULL, NULL, shaderFlags, 0, &compiledFX, &errors);
    if( FAILED(hr) )
    {
        char *errMsg = static_cast<char *>(errors->GetBufferPointer());
        errMsg[errors->GetBufferSize()-1] = '\0';
        MessageBoxA( NULL, errMsg, "Effect compilation failed", MB_OK|MB_ICONERROR );
        errors->Release();
        return hr;
    }
    hr = D3D10CreateEffectFromMemory(compiledFX->GetBufferPointer(), compiledFX->GetBufferSize(), 
                                     0, g_D3DDevice, NULL, &g_Effect);
    compiledFX->Release();

    if( FAILED( hr ) )
    {
        MessageBox( NULL, L"Effect creation failed", L"Error", MB_OK );
        return hr;
    }

    // Obtain the technique
    g_Technique = g_Effect->GetTechniqueByName("Render");

    // Define the input layout
    D3D10_INPUT_ELEMENT_DESC layout[] =
    {
        { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D10_INPUT_PER_VERTEX_DATA, 0 },  
        { "COLOR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 3*sizeof(float), D3D10_INPUT_PER_VERTEX_DATA, 0 }
    };

    // Create the input layout
    D3D10_PASS_DESC passDesc;
    g_Technique->GetPassByIndex(0)->GetDesc(&passDesc);
    hr = g_D3DDevice->CreateInputLayout(layout, sizeof(layout)/sizeof(layout[0]), passDesc.pIAInputSignature,
                                        passDesc.IAInputSignatureSize, &g_VertexLayout);
    if( FAILED( hr ) )
        return hr;
    const int VERTEX_SIZE = (3 + 4)*sizeof(float); // 3 floats for POSITION + 4 floats for COLOR

    // Set the input layout
    g_D3DDevice->IASetInputLayout(g_VertexLayout);

    // Create vertex buffer
    D3D10_BUFFER_DESC bd;
    bd.Usage = D3D10_USAGE_DYNAMIC;
    bd.ByteWidth = VERTEX_SIZE * NB_VERTS;
    bd.BindFlags = D3D10_BIND_VERTEX_BUFFER;
    bd.CPUAccessFlags = D3D10_CPU_ACCESS_WRITE;
    bd.MiscFlags = 0;
    hr = g_D3DDevice->CreateBuffer(&bd, NULL, &g_VertexBuffer);
    if( FAILED( hr ) )
        return hr;

    // Set vertex buffer
    UINT stride = VERTEX_SIZE;
    UINT offset = 0;
    g_D3DDevice->IASetVertexBuffers(0, 1, &g_VertexBuffer, &stride, &offset);

    // Set primitive topology
    g_D3DDevice->IASetPrimitiveTopology(D3D10_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

    // Blend state
    D3D10_BLEND_DESC bsd;
    bsd.AlphaToCoverageEnable = FALSE;
    for(int i=0; i<8; ++i)
    {
        bsd.BlendEnable[i] = TRUE;
        bsd.RenderTargetWriteMask[i] = D3D10_COLOR_WRITE_ENABLE_ALL;
    }
    bsd.SrcBlend = D3D10_BLEND_SRC_ALPHA;
    bsd.DestBlend = D3D10_BLEND_INV_SRC_ALPHA;
    bsd.BlendOp =  D3D10_BLEND_OP_ADD;
    bsd.SrcBlendAlpha = D3D10_BLEND_SRC_ALPHA;
    bsd.DestBlendAlpha = D3D10_BLEND_INV_SRC_ALPHA;
    bsd.BlendOpAlpha = D3D10_BLEND_OP_ADD;
    g_D3DDevice->CreateBlendState(&bsd, &g_BlendState);
    float blendFactors[4] = { 1, 1, 1, 1 };
    g_D3DDevice->OMSetBlendState(g_BlendState, blendFactors, 0xffffffff);

    // Rasterizer state
    D3D10_RASTERIZER_DESC rs;
    ZeroMemory(&rs, sizeof(rs));
    rs.FillMode = D3D10_FILL_SOLID;
    rs.CullMode = D3D10_CULL_NONE;
    g_D3DDevice->CreateRasterizerState(&rs, &g_RasterState);
    g_D3DDevice->RSSetState(g_RasterState);

    return S_OK;
}


// Clean up D3D objects
void Cleanup()
{
    if( g_D3DDevice ) 
        g_D3DDevice->ClearState();

    if( g_BlendState )
        g_BlendState->Release();
    if( g_RasterState )
        g_RasterState->Release();
    if( g_VertexBuffer ) 
        g_VertexBuffer->Release();
    if( g_VertexLayout ) 
        g_VertexLayout->Release();
    if( g_Effect ) 
        g_Effect->Release();
    if( g_RenderTargetView ) 
        g_RenderTargetView->Release();
    if( g_SwapChain && g_SwapChainDesc.Windowed ) 
        g_SwapChain->Release();
    if( g_D3DDevice ) 
        g_D3DDevice->Release();
}


// Called every time the application receives a message
LRESULT CALLBACK MessageProc(HWND wnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    // Send event message to AntTweakBar
    if( TwEventWin(wnd, message, wParam, lParam) ) 
    {
        // Event has been handled by AntTweakBar.
        // Note: TwEventWin send a WM_PAINT in this case and
        //       the window will be repaint next time.
        return 0;
    }

    switch(message) 
    {
        case WM_PAINT:
            {
                PAINTSTRUCT ps;
                BeginPaint(wnd, &ps);

                if( g_D3DDevice )
                    Render();

                EndPaint(wnd, &ps);
            }
            return 0;
        case WM_SIZE:   // Window size has been changed
            if( g_D3DDevice )   // Resize D3D render target
            {
                // Release render target
                ID3D10RenderTargetView *nullRTV = NULL;
                g_D3DDevice->OMSetRenderTargets(1, &nullRTV, NULL);
                if( g_RenderTargetView )
                {
                    g_RenderTargetView->Release();
                    g_RenderTargetView = NULL;
                }

                if( g_SwapChain )
                {
                    // Resize swap chain
                    g_SwapChainDesc.BufferDesc.Width = LOWORD(lParam);
                    g_SwapChainDesc.BufferDesc.Height = HIWORD(lParam);
                    g_SwapChain->ResizeBuffers(g_SwapChainDesc.BufferCount, g_SwapChainDesc.BufferDesc.Width, 
                        g_SwapChainDesc.BufferDesc.Height, g_SwapChainDesc.BufferDesc.Format, g_SwapChainDesc.Flags);

                    // Re-create a render target view
                    ID3D10Texture2D *buffer;
                    g_SwapChain->GetBuffer(0, __uuidof( ID3D10Texture2D ), (LPVOID*)&buffer);
                    g_D3DDevice->CreateRenderTargetView(buffer, NULL, &g_RenderTargetView);
                    buffer->Release();
                    g_D3DDevice->OMSetRenderTargets(1, &g_RenderTargetView, NULL);

                    // Setup the viewport
                    D3D10_VIEWPORT vp;
                    vp.Width = g_SwapChainDesc.BufferDesc.Width;
                    vp.Height = g_SwapChainDesc.BufferDesc.Height;
                    vp.MinDepth = 0.0f;
                    vp.MaxDepth = 1.0f;
                    vp.TopLeftX = 0;
                    vp.TopLeftY = 0;
                    g_D3DDevice->RSSetViewports(1, &vp);
                }

                // TwWindowSize has been called by TwEventWin32, so it is not necessary to call it again here.
            }
            return 0;
        case WM_CHAR:
            if( wParam==VK_ESCAPE )
                PostQuitMessage(0);
            return 0;
        case WM_DESTROY:
            PostQuitMessage(0);
            return 0;
        default:
            return DefWindowProc(wnd, message, wParam, lParam);
    }
}


// Render scene
void Render()
{
    // Update the vertex buffer
    float a = (float)g_Angle * (3.14159265358979f/180.0f);
    float ca = cosf(a);
    float sa = sinf(a);
    float ratio = (float)g_SwapChainDesc.BufferDesc.Height/g_SwapChainDesc.BufferDesc.Width;
    float *vertices;
    g_VertexBuffer->Map(D3D10_MAP_WRITE_DISCARD, 0, (void **)&vertices);
    // Vertex 0 (3 floats for POSITION + 4 floats for COLOR)
    vertices[0*7 + 0] = g_Scale*(ca*g_Positions[0].X - sa*g_Positions[0].Y)*ratio;
    vertices[0*7 + 1] = g_Scale*(sa*g_Positions[0].X + ca*g_Positions[0].Y);
    vertices[0*7 + 2] = 0;
    vertices[0*7 + 3] = g_Colors[0][0];
    vertices[0*7 + 4] = g_Colors[0][1];
    vertices[0*7 + 5] = g_Colors[0][2];
    vertices[0*7 + 6] = g_Colors[0][3];
    // Vertex 1 (3 floats for POSITION + 4 floats for COLOR)
    vertices[1*7 + 0] = g_Scale*(ca*g_Positions[1].X - sa*g_Positions[1].Y)*ratio;
    vertices[1*7 + 1] = g_Scale*(sa*g_Positions[1].X + ca*g_Positions[1].Y);
    vertices[1*7 + 2] = 0;
    vertices[1*7 + 3] = g_Colors[1][0];
    vertices[1*7 + 4] = g_Colors[1][1];
    vertices[1*7 + 5] = g_Colors[1][2];
    vertices[1*7 + 6] = g_Colors[1][3];
    // Vertex 2 (3 floats for POSITION + 4 floats for COLOR)
    vertices[2*7 + 0] = g_Scale*(ca*g_Positions[2].X - sa*g_Positions[2].Y)*ratio;
    vertices[2*7 + 1] = g_Scale*(sa*g_Positions[2].X + ca*g_Positions[2].Y);
    vertices[2*7 + 2] = 0;
    vertices[2*7 + 3] = g_Colors[2][0];
    vertices[2*7 + 4] = g_Colors[2][1];
    vertices[2*7 + 5] = g_Colors[2][2];
    vertices[2*7 + 6] = g_Colors[2][3];
    g_VertexBuffer->Unmap();

    // Clear the back buffer 
    float clearColor[4] = {0.125f, 0.125f, 0.3f, 1.0f};
    g_D3DDevice->ClearRenderTargetView(g_RenderTargetView, clearColor);

    // Render a triangle
    D3D10_TECHNIQUE_DESC techDesc;
    g_Technique->GetDesc(&techDesc);
    for( UINT p=0; p<techDesc.Passes; ++p )
    {
        g_Technique->GetPassByIndex(p)->Apply(0);
        g_D3DDevice->Draw(3, 0);
    }

	// Draw tweak bars
    TwDraw();

    // Present the information rendered to the back buffer to the front buffer (the screen)
    g_SwapChain->Present(0, 0);
}


