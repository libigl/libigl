//  ---------------------------------------------------------------------------
//
//  @file       MiniSFML16.h
//  @brief      A subset of SFML 1.6 definitions needed to compile helper
//              functions implemented in TwEventSFML.cpp
//
//  notes:    - Private header
//            - AntTweakBar.dll does not need to link with SFML, 
//              it just needs some definitions for its helper functions.
//            - This header is provided to avoid the need of having SFML
//              installed to recompile AntTweakBar.
//              It declares a small and incomplete part of SFML classes.
//              For instance, many non-virtual methods have been stripped out.
//            - Do not use this header in your own programs, better use the
//              SFML headers from the actual SFML library SDK :
//              http://www.sfml-dev.org
//
//  ---------------------------------------------------------------------------

#if !defined MINI_SFML16_INCLUDED
#define MINI_SFML16_INCLUDED

namespace sf 
{

    namespace Key { enum Code 
    {
        A = 'a', B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z,
        Num0 = '0', Num1, Num2, Num3, Num4, Num5, Num6, Num7, Num8, Num9,
        Escape = 256,
        LControl, LShift, LAlt, LSystem,
        RControl, RShift, RAlt, RSystem, Menu,
        LBracket, RBracket, SemiColon, Comma, Period, Quote, Slash, BackSlash,
        Tilde, Equal, Dash, Space, Return, Back, Tab,
        PageUp, PageDown, End, Home, Insert, Delete,
        Add, Subtract, Multiply, Divide, Left, Right, Up, Down,
        Numpad0, Numpad1, Numpad2, Numpad3, Numpad4, Numpad5, Numpad6, Numpad7, Numpad8, Numpad9,
        F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, Pause,
        Count
    }; }

    namespace Mouse { enum Button { Left, Right, Middle, XButton1, XButton2, Count }; }

    namespace Joy { enum Axis { AxisX, AxisY, AxisZ, AxisR, AxisU, AxisV, AxisPOV, AxisCount }; 
                    enum { Count = 4, ButtonCount = 32}; }

    typedef unsigned char Uint8;
    typedef unsigned int Uint32;

    class Event
    {
    public :
        struct KeyEvent { Key::Code Code; bool Alt, Control, Shift; };
        struct TextEvent { Uint32 Unicode; };
        struct MouseMoveEvent { int X, Y; };
        struct MouseButtonEvent { Mouse::Button Button; int X, Y; };
        struct MouseWheelEvent { int Delta; };
        struct JoyMoveEvent { unsigned int JoystickId; Joy::Axis Axis; float Position; };
        struct JoyButtonEvent { unsigned int JoystickId, Button; };
        struct SizeEvent { unsigned int Width, Height; };
        enum EventType
        {
            Closed, Resized, LostFocus, GainedFocus, TextEntered, KeyPressed, KeyReleased,
            MouseWheelMoved, MouseButtonPressed, MouseButtonReleased, MouseMoved, MouseEntered, MouseLeft,
            JoyButtonPressed, JoyButtonReleased, JoyMoved
        };
        EventType Type;
        union
        {
            KeyEvent         Key;
            TextEvent        Text;
            MouseMoveEvent   MouseMove;
            MouseButtonEvent MouseButton;
            MouseWheelEvent  MouseWheel;
            JoyMoveEvent     JoyMove;
            JoyButtonEvent   JoyButton;
            SizeEvent        Size;
        };
    };

} // namespace sf


#ifdef USE_MINI_SFML 
// we also need some the definition of sf::RenderWindow to compile our SFML example

#include <string>
#include <queue>

namespace sf 
{

    class Input;
    class Drawable;
    typedef void* WindowHandle;
    namespace priv { class WindowImpl; }

    class WindowListener
    {
    public :
        virtual void OnEvent(const Event& EventReceived) = 0;
    protected :
        virtual ~WindowListener();
    };

    class VideoMode
    {
    public :
        VideoMode(unsigned int ModeWidth, unsigned int ModeHeight, unsigned int ModeBpp = 32);
        unsigned int Width, Height, BitsPerPixel;
    };

    namespace Style { enum { None = 0, Titlebar = 1 << 0, Resize = 1 << 1, Close = 1 << 2, Fullscreen = 1 << 3 }; }

    struct WindowSettings
    {
        explicit WindowSettings(unsigned int Depth = 24, unsigned int Stencil = 8, unsigned int Antialiasing = 0);
        unsigned int DepthBits, StencilBits, AntialiasingLevel;
    };

    class Clock
    {
    public :
        Clock();
        float GetElapsedTime() const;
        void Reset();
    private :
        double myStartTime;
    };

    class Input : public WindowListener
    {
    private :
        virtual void OnEvent(const Event& EventReceived);
        bool  myKeys[Key::Count];           
        bool  myMouseButtons[Mouse::Count]; 
        int   myMouseX;                     
        int   myMouseY;                     
        bool  myJoystickButtons[Joy::Count][Joy::ButtonCount];     
        float myJoystickAxis[Joy::Count][Joy::AxisCount];
    };

    class Window : public WindowListener
    {
    public :
        Window(VideoMode Mode, const std::string& Title, unsigned long WindowStyle = Style::Resize | Style::Close, const WindowSettings& Params = WindowSettings());
        virtual ~Window();
        void Close();
        bool IsOpened() const;
        unsigned int GetWidth() const;
        unsigned int GetHeight() const;
        bool GetEvent(Event& EventReceived);
        void Display();
    private :
        virtual void OnCreate();
        virtual void OnEvent(const Event& EventReceived);
        priv::WindowImpl* myWindow;        
        std::queue<Event> myEvents;        
        Input             myInput;         
        Clock             myClock;         
        WindowSettings    mySettings;      
        float             myLastFrameTime; 
        bool              myIsExternal;    
        unsigned int      myFramerateLimit;
        int               mySetCursorPosX; 
        int               mySetCursorPosY; 
    };

    template <typename T> class Vector2 { public : T x, y; };
    typedef Vector2<float> Vector2f;

    template <typename T> class Rect { public : T Left, Top, Right, Bottom; };
    typedef Rect<float> FloatRect;

    class Matrix3 { private : float myData[16]; };

    class View
    {
    private :
        sf::Vector2f myCenter;
        sf::Vector2f myHalfSize;
        FloatRect    myRect;
        Matrix3      myMatrix;
        bool         myNeedUpdate;
    };

    class RenderTarget
    {
    public :
        virtual ~RenderTarget();
        virtual void Draw(const Drawable& Object);
        virtual unsigned int GetWidth() const = 0;
        virtual unsigned int GetHeight() const = 0;
        void PreserveOpenGLStates(bool Preserve);
    private :
        virtual bool Activate(bool Active) = 0;
        View        myDefaultView;
        const View* myCurrentView;
        bool        myPreserveStates;
        bool        myIsDrawing;
    };

    class RenderWindow : public Window, public RenderTarget
    {
    public :
        RenderWindow(VideoMode Mode, const std::string& Title, unsigned long WindowStyle = Style::Resize | Style::Close, const WindowSettings& Params = WindowSettings());
        virtual ~RenderWindow();
        virtual unsigned int GetWidth() const;
        virtual unsigned int GetHeight() const;
    private :
        virtual void OnCreate();
        virtual bool Activate(bool Active);
    };

} // namespace sf

#endif // USE_MINI_SFML

#endif // !defined MINI_SFML16_INCLUDED

