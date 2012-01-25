//  ---------------------------------------------------------------------------
//
//  @file       MiniSDL13.h
//  @brief      A subset of SDL 1.3 definitions needed to compile helper 
//              functions implemented in TwEventSDL13.c
//
//  notes:    - Private header
//            - AntTweakBar.dll does not need to link with SDL, 
//              it just needs some definitions for its helper functions.
//            - This header is provided to avoid the need of having SDL
//              installed to recompile AntTweakBar.
//            - Do not use this header in your own programs, better use the
//              SDL.h header from the actual SDL library SDK :
//              http://www.libsdl.org
//
//  ---------------------------------------------------------------------------

#if !defined MINI_SDL_INCLUDED
#define MINI_SDL_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#define SDL_MAJOR_VERSION	1
#define SDL_MINOR_VERSION	3

#if defined(_WIN32) || defined(_WIN64)
#   define SDL_DECLSPEC __declspec(dllimport)
#   define SDL_CALL __cdecl
#else
#   define SDL_DECLSPEC
#   define SDL_CALL
#endif

typedef unsigned char   Uint8;
typedef signed char     Sint8;
typedef unsigned short  Uint16;
typedef signed short    Sint16;
typedef unsigned int    Uint32;
typedef signed int      Sint32;

#define SDLK_SCANCODE_MASK (1<<30)
#define SDL_SCANCODE_TO_KEYCODE(X)	(X | SDLK_SCANCODE_MASK)

// Subset of SDL scancodes.
// Note: some SDL scancodes seems to be wrong in the original
// SDL scancode header file. 
typedef enum {
    SDL_SCANCODE_F1         = 58,
    SDL_SCANCODE_F2         = 59,
    SDL_SCANCODE_F3         = 60,
    SDL_SCANCODE_F4         = 61,
    SDL_SCANCODE_F5         = 62,
    SDL_SCANCODE_F6         = 63,
    SDL_SCANCODE_F7         = 64,
    SDL_SCANCODE_F8         = 65,
    SDL_SCANCODE_F9         = 66,
    SDL_SCANCODE_F10        = 67,
    SDL_SCANCODE_F11        = 68,
    SDL_SCANCODE_F12        = 69,
    SDL_SCANCODE_INSERT     = 98, //73,
    SDL_SCANCODE_HOME       = 95, //74,
    SDL_SCANCODE_PAGEUP     = 97, //75,
    SDL_SCANCODE_DELETE     = 99, //76,
    SDL_SCANCODE_END        = 89, //77,
    SDL_SCANCODE_PAGEDOWN   = 91, //78,
    SDL_SCANCODE_RIGHT      = 94, //79,
    SDL_SCANCODE_LEFT       = 92, //80,
    SDL_SCANCODE_DOWN       = 90, //81,
    SDL_SCANCODE_UP         = 96  //82
} SDL_scancode;

// Subset of SDL keysym
typedef enum {
    SDLK_BACKSPACE  = 8,
    SDLK_TAB        = 9,
    SDLK_CLEAR      = 12,
    SDLK_RETURN     = 13,
    SDLK_PAUSE      = 19,
    SDLK_ESCAPE     = 27,
    SDLK_DELETE     = 127,
    SDLK_UP         = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_UP),
    SDLK_DOWN       = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_DOWN),
    SDLK_RIGHT      = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_RIGHT),
    SDLK_LEFT       = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_LEFT),
    SDLK_INSERT     = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_INSERT),
    SDLK_HOME       = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_HOME),
    SDLK_END        = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_END),
    SDLK_PAGEUP     = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_PAGEUP),
    SDLK_PAGEDOWN   = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_PAGEDOWN),
    SDLK_F1         = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_F1),
    SDLK_F2         = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_F2),
    SDLK_F3         = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_F3),
    SDLK_F4         = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_F4),
    SDLK_F5         = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_F5),
    SDLK_F6         = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_F6),
    SDLK_F7         = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_F7),
    SDLK_F8         = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_F8),
    SDLK_F9         = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_F9),
    SDLK_F10        = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_F10),
    SDLK_F11        = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_F11),
    SDLK_F12        = SDL_SCANCODE_TO_KEYCODE(SDL_SCANCODE_F12)
} SDLKey;

typedef enum {
    KMOD_NONE       = 0x0000,
    KMOD_LSHIFT     = 0x0001,
    KMOD_RSHIFT     = 0x0002,
    KMOD_LCTRL      = 0x0040,
    KMOD_RCTRL      = 0x0080,
    KMOD_LALT       = 0x0100,
    KMOD_RALT       = 0x0200,
    KMOD_LGUI       = 0x0400,
    KMOD_RGUI       = 0x0800,
    KMOD_NUM        = 0x1000,
    KMOD_CAPS       = 0x2000,
    KMOD_MODE       = 0x4000,
    KMOD_RESERVED   = 0x8000
} SDLMod;

#define KMOD_CTRL   (KMOD_LCTRL|KMOD_RCTRL)
#define KMOD_SHIFT  (KMOD_LSHIFT|KMOD_RSHIFT)
#define KMOD_ALT    (KMOD_LALT|KMOD_RALT)
#define KMOD_GUI    (KMOD_LGUI|KMOD_RGUI)

typedef enum { 
    SDL_NOEVENT     = 0,
    SDL_WINDOWEVENT,    
    SDL_KEYDOWN,        
    SDL_KEYUP,          
    SDL_TEXTEDITING,    
    SDL_TEXTINPUT,      
    SDL_MOUSEMOTION,    
    SDL_MOUSEBUTTONDOWN,
    SDL_MOUSEBUTTONUP,  
    SDL_MOUSEWHEEL,     
    SDL_JOYAXISMOTION,  
    SDL_JOYBALLMOTION,  
    SDL_JOYHATMOTION,   
    SDL_JOYBUTTONDOWN,  
    SDL_JOYBUTTONUP,    
    SDL_QUIT,           
    SDL_SYSWMEVENT,     
    SDL_PROXIMITYIN,    
    SDL_PROXIMITYOUT,   
    SDL_EVENT_RESERVED1,
    SDL_EVENT_RESERVED2,
    SDL_EVENT_RESERVED3,
    SDL_USEREVENT = 24,
    SDL_NUMEVENTS = 32
} SDL_EventType;

#define SDL_ACTIVEEVENT	SDL_EVENT_RESERVED1
#define SDL_VIDEORESIZE	SDL_EVENT_RESERVED2
#define SDL_VIDEOEXPOSE	SDL_EVENT_RESERVED3

typedef Uint32 SDL_WindowID;

typedef struct SDL_keysym {
    SDL_scancode scancode;
    SDLKey sym;
    Uint16 mod;
    Uint32 unicode;
} SDL_keysym;

typedef struct SDL_WindowEvent {
    Uint8 type;             
    SDL_WindowID windowID;  
    Uint8 event;            
    int data1;              
    int data2;              
} SDL_WindowEvent;

typedef struct SDL_KeyboardEvent {
    Uint8 type;            
    SDL_WindowID windowID; 
    Uint8 which;           
    Uint8 state;           
    SDL_keysym keysym;     
} SDL_KeyboardEvent;

#define SDL_TEXTEDITINGEVENT_TEXT_SIZE (32)
typedef struct SDL_TextEditingEvent {
    Uint8 type;                                
    char text[SDL_TEXTEDITINGEVENT_TEXT_SIZE]; 
    int start;                                 
    int length;                                
} SDL_TextEditingEvent;

#define SDL_TEXTINPUTEVENT_TEXT_SIZE (32)
typedef struct SDL_TextInputEvent {
    Uint8 type;                              
    SDL_WindowID windowID;                   
    Uint8 which;                             
    char text[SDL_TEXTINPUTEVENT_TEXT_SIZE];
} SDL_TextInputEvent;

typedef struct SDL_MouseMotionEvent {
    Uint8 type;            
    SDL_WindowID windowID; 
    Uint8 which;           
    Uint8 state;           
    int x;                 
    int y;                 
    int z;                 
    int pressure;          
    int pressure_max;      
    int pressure_min;      
    int rotation;          
    int tilt;              
    int cursor;            
    int xrel;              
    int yrel;              
} SDL_MouseMotionEvent;

typedef struct SDL_MouseButtonEvent {
    Uint8 type;            
    SDL_WindowID windowID; 
    Uint8 which;           
    Uint8 button;          
    Uint8 state;           
    int x;                 
    int y;                 
} SDL_MouseButtonEvent;

typedef struct SDL_MouseWheelEvent {
    Uint8 type;            
    SDL_WindowID windowID; 
    Uint8 which;           
    int x;                 
    int y;                 
} SDL_MouseWheelEvent;

typedef struct SDL_JoyAxisEvent {
    Uint8 type;
    Uint8 which;
    Uint8 axis;
    Sint16 value;
} SDL_JoyAxisEvent;

typedef struct SDL_JoyBallEvent {
    Uint8 type;
    Uint8 which;
    Uint8 ball;
    Sint16 xrel;
    Sint16 yrel;
} SDL_JoyBallEvent;

typedef struct SDL_JoyHatEvent {
    Uint8 type;
    Uint8 which;
    Uint8 hat;
    Uint8 value;
} SDL_JoyHatEvent;

typedef struct SDL_JoyButtonEvent {
    Uint8 type;
    Uint8 which;
    Uint8 button;
    Uint8 state;
} SDL_JoyButtonEvent;

typedef struct SDL_ActiveEvent
{
    Uint8 type;
    Uint8 gain;
    Uint8 state;
} SDL_ActiveEvent;

typedef struct SDL_ResizeEvent {
    Uint8 type;
    int w;
    int h;
} SDL_ResizeEvent;

typedef struct SDL_ExposeEvent {
    Uint8 type;
} SDL_ExposeEvent;

typedef struct SDL_QuitEvent {
    Uint8 type;
} SDL_QuitEvent;

typedef struct SDL_UserEvent {
    Uint8 type;
    SDL_WindowID windowID; 
    int code;
    void *data1;
    void *data2;
} SDL_UserEvent;

struct SDL_SysWMmsg;
typedef struct SDL_SysWMmsg SDL_SysWMmsg;
typedef struct SDL_SysWMEvent {
    Uint8 type;
    SDL_SysWMmsg *msg;
} SDL_SysWMEvent;

typedef struct SDL_ProximityEvent
{
    Uint8 type;
    SDL_WindowID windowID;
    Uint8 which;
    int cursor;
    int x;
    int y;
} SDL_ProximityEvent;

typedef union SDL_Event {
    Uint8 type;                    
    SDL_WindowEvent window;        
    SDL_KeyboardEvent key;         
    SDL_TextEditingEvent edit;     
    SDL_TextInputEvent text;       
    SDL_MouseMotionEvent motion;   
    SDL_MouseButtonEvent button;   
    SDL_MouseWheelEvent wheel;     
    SDL_JoyAxisEvent jaxis;        
    SDL_JoyBallEvent jball;        
    SDL_JoyHatEvent jhat;          
    SDL_JoyButtonEvent jbutton;    
    SDL_QuitEvent quit;            
    SDL_UserEvent user;            
    SDL_SysWMEvent syswm;          
    SDL_ProximityEvent proximity;  
    SDL_ActiveEvent active;
    SDL_ResizeEvent resize;
} SDL_Event;

typedef struct SDL_PixelFormat {
    void  *palette;
    Uint8  BitsPerPixel;
    Uint8  BytesPerPixel;
    Uint8  Rloss;
    Uint8  Gloss;
    Uint8  Bloss;
    Uint8  Aloss;
    Uint8  Rshift;
    Uint8  Gshift;
    Uint8  Bshift;
    Uint8  Ashift;
    Uint32 Rmask;
    Uint32 Gmask;
    Uint32 Bmask;
    Uint32 Amask;
} SDL_PixelFormat;

typedef enum SDL_GLattr {
    SDL_GL_RED_SIZE,
    SDL_GL_GREEN_SIZE,
    SDL_GL_BLUE_SIZE,
    SDL_GL_ALPHA_SIZE,
    SDL_GL_BUFFER_SIZE,
    SDL_GL_DOUBLEBUFFER,
    SDL_GL_DEPTH_SIZE,
    SDL_GL_STENCIL_SIZE,
    SDL_GL_ACCUM_RED_SIZE,
    SDL_GL_ACCUM_GREEN_SIZE,
    SDL_GL_ACCUM_BLUE_SIZE,
    SDL_GL_ACCUM_ALPHA_SIZE,
    SDL_GL_STEREO,
    SDL_GL_MULTISAMPLEBUFFERS,
    SDL_GL_MULTISAMPLESAMPLES,
    SDL_GL_ACCELERATED_VISUAL,
    SDL_GL_RETAINED_BACKING,
    SDL_GL_CONTEXT_MAJOR_VERSION,
    SDL_GL_CONTEXT_MINOR_VERSION
} SDL_GLattr;

typedef struct SDL_VideoInfo {
    Uint32 hw_available :1;
    Uint32 wm_available :1;
    Uint32 UnusedBits1  :6;
    Uint32 UnusedBits2  :1;
    Uint32 blit_hw      :1;
    Uint32 blit_hw_CC   :1;
    Uint32 blit_hw_A    :1;
    Uint32 blit_sw      :1;
    Uint32 blit_sw_CC   :1;
    Uint32 blit_sw_A    :1;
    Uint32 blit_fill    :1;
    Uint32 UnusedBits3  :16;
    Uint32 video_mem;
    SDL_PixelFormat *vfmt;
    int current_w;
    int current_h;
} SDL_VideoInfo;

#define SDL_INIT_VIDEO  0x00000020

#define SDL_ANYFORMAT   0x00100000
#define SDL_HWPALETTE   0x00200000
#define SDL_DOUBLEBUF   0x00400000
#define SDL_FULLSCREEN  0x00800000
#define SDL_RESIZABLE   0x01000000
#define SDL_NOFRAME     0x02000000
#define SDL_OPENGL      0x04000000
#define SDL_HWSURFACE   0x08000001

#define SDL_DEFAULT_REPEAT_DELAY    500
#define SDL_DEFAULT_REPEAT_INTERVAL 30


// functions subset
extern SDL_DECLSPEC int    SDL_CALL SDL_Init(Uint32 flags);
extern SDL_DECLSPEC void   SDL_CALL SDL_Quit();
extern SDL_DECLSPEC char * SDL_CALL SDL_GetError();
extern SDL_DECLSPEC const SDL_VideoInfo * SDL_CALL SDL_GetVideoInfo();
extern SDL_DECLSPEC struct SDL_Surface * SDL_CALL SDL_SetVideoMode(int width, int height, int bpp, Uint32 flags);
extern SDL_DECLSPEC int    SDL_CALL SDL_GL_SetAttribute(SDL_GLattr attr, int value);
extern SDL_DECLSPEC void   SDL_CALL SDL_GL_SwapBuffers();
extern SDL_DECLSPEC void   SDL_CALL SDL_WM_SetCaption(const char *title, const char *icon);
extern SDL_DECLSPEC void   SDL_CALL SDL_WM_GetCaption(char **title, char **icon);
extern SDL_DECLSPEC int    SDL_CALL SDL_EnableUNICODE(int enable);
extern SDL_DECLSPEC int    SDL_CALL SDL_EnableKeyRepeat(int delay, int interval);
extern SDL_DECLSPEC Uint32 SDL_CALL SDL_GetTicks();
extern SDL_DECLSPEC int    SDL_CALL SDL_PollEvent(SDL_Event *event);
extern SDL_DECLSPEC int    SDL_CALL SDL_WaitEvent(SDL_Event *event);
extern SDL_DECLSPEC int    SDL_CALL SDL_PushEvent(SDL_Event *event);


#ifdef __cplusplus
}
#endif

#endif // !defined MINI_SDL_INCLUDED

