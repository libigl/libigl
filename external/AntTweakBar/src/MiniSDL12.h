//  ---------------------------------------------------------------------------
//
//  @file       MiniSDL12.h
//  @brief      A subset of SDL 1.2 definitions needed to compile helper
//              functions implemented in TwEventSDL12.c
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

#if !defined MINI_SDL12_INCLUDED
#define MINI_SDL12_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#define SDL_MAJOR_VERSION	1
#define SDL_MINOR_VERSION	2

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

// Subset of SDL keysym
typedef enum {
    SDLK_BACKSPACE  = 8,
    SDLK_TAB        = 9,
    SDLK_CLEAR      = 12,
    SDLK_RETURN     = 13,
    SDLK_PAUSE      = 19,
    SDLK_ESCAPE     = 27,
    SDLK_DELETE     = 127,
    SDLK_UP         = 273,
    SDLK_DOWN       = 274,
    SDLK_RIGHT      = 275,
    SDLK_LEFT       = 276,
    SDLK_INSERT     = 277,
    SDLK_HOME       = 278,
    SDLK_END        = 279,
    SDLK_PAGEUP     = 280,
    SDLK_PAGEDOWN   = 281,
    SDLK_F1         = 282,
    SDLK_F2         = 283,
    SDLK_F3         = 284,
    SDLK_F4         = 285,
    SDLK_F5         = 286,
    SDLK_F6         = 287,
    SDLK_F7         = 288,
    SDLK_F8         = 289,
    SDLK_F9         = 290,
    SDLK_F10        = 291,
    SDLK_F11        = 292,
    SDLK_F12        = 293,
} SDLKey;

typedef enum {
    KMOD_NONE       = 0x0000,
    KMOD_LSHIFT     = 0x0001,
    KMOD_RSHIFT     = 0x0002,
    KMOD_LCTRL      = 0x0040,
    KMOD_RCTRL      = 0x0080,
    KMOD_LALT       = 0x0100,
    KMOD_RALT       = 0x0200,
    KMOD_LMETA      = 0x0400,
    KMOD_RMETA      = 0x0800,
    KMOD_NUM        = 0x1000,
    KMOD_CAPS       = 0x2000,
    KMOD_MODE       = 0x4000,
    KMOD_RESERVED   = 0x8000
} SDLMod;

#define KMOD_CTRL   (KMOD_LCTRL|KMOD_RCTRL)
#define KMOD_SHIFT  (KMOD_LSHIFT|KMOD_RSHIFT)
#define KMOD_ALT    (KMOD_LALT|KMOD_RALT)
#define KMOD_META   (KMOD_LMETA|KMOD_RMETA)

typedef enum { 
    SDL_NOEVENT = 0,
    SDL_ACTIVEEVENT,
    SDL_KEYDOWN,
    SDL_KEYUP,
    SDL_MOUSEMOTION,
    SDL_MOUSEBUTTONDOWN,
    SDL_MOUSEBUTTONUP,
    SDL_JOYAXISMOTION,
    SDL_JOYBALLMOTION,
    SDL_JOYHATMOTION,
    SDL_JOYBUTTONDOWN,
    SDL_JOYBUTTONUP,
    SDL_QUIT,
    SDL_SYSWMEVENT,
    SDL_EVENT_RESERVEDA,
    SDL_EVENT_RESERVEDB,
    SDL_VIDEORESIZE,
    SDL_VIDEOEXPOSE,
    SDL_EVENT_RESERVED2,
    SDL_EVENT_RESERVED3,
    SDL_EVENT_RESERVED4,
    SDL_EVENT_RESERVED5,
    SDL_EVENT_RESERVED6,
    SDL_EVENT_RESERVED7,
    SDL_USEREVENT = 24,
    SDL_NUMEVENTS = 32
} SDLEventEnum;

typedef struct SDL_keysym {
    Uint8 scancode;
    SDLKey sym;
    SDLMod mod;
    Uint16 unicode;
} SDL_keysym;

typedef struct SDL_ActiveEvent {
    Uint8 type;
    Uint8 gain;
    Uint8 state;
} SDL_ActiveEvent;

typedef struct SDL_KeyboardEvent {
    Uint8 type;
    Uint8 which;
    Uint8 state;
    SDL_keysym keysym;
} SDL_KeyboardEvent;

typedef struct SDL_MouseMotionEvent {
    Uint8 type;
    Uint8 which;
    Uint8 state;
    Uint16 x, y;
    Sint16 xrel;
    Sint16 yrel;
} SDL_MouseMotionEvent;

typedef struct SDL_MouseButtonEvent {
    Uint8 type;
    Uint8 which;
    Uint8 button;
    Uint8 state;
    Uint16 x, y;
} SDL_MouseButtonEvent;

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

typedef union {
    Uint8 type;
    SDL_ActiveEvent active;
    SDL_KeyboardEvent key;
    SDL_MouseMotionEvent motion;
    SDL_MouseButtonEvent button;
    SDL_JoyAxisEvent jaxis;
    SDL_JoyBallEvent jball;
    SDL_JoyHatEvent jhat;
    SDL_JoyButtonEvent jbutton;
    SDL_ResizeEvent resize;
    SDL_ExposeEvent expose;
    SDL_QuitEvent quit;
    SDL_UserEvent user;
    SDL_SysWMEvent syswm;
char full[56];
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
    Uint32 colorkey;
    Uint8  alpha;
} SDL_PixelFormat;

typedef enum {
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

#define SDL_SWSURFACE   0x00000000
#define SDL_HWSURFACE   0x00000001
#define SDL_ASYNCBLIT   0x00000004
#define SDL_ANYFORMAT   0x10000000
#define SDL_HWPALETTE   0x20000000
#define SDL_DOUBLEBUF   0x40000000
#define SDL_FULLSCREEN  0x80000000
#define SDL_OPENGL      0x00000002
#define SDL_OPENGLBLIT  0x0000000A
#define SDL_RESIZABLE   0x00000010
#define SDL_NOFRAME     0x00000020

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

#endif // !defined MINI_SDL12_INCLUDED

