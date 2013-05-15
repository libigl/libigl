//  ---------------------------------------------------------------------------
//
//  @file       TwEventSFML.cpp
//  @brief      Helper: 
//              translate and re-send mouse and keyboard events 
//              from SFML 1.6 event loop to AntTweakBar
//  
//  @author     Philippe Decaudin
//  @license    This file is part of the AntTweakBar library.
//              For conditions of distribution and use, see License.txt
//
//  ---------------------------------------------------------------------------


#include "MiniSFML16.h" // a subset of SFML 1.6 headers needed to compile TwEventSFML.cpp
// note: AntTweakBar.dll does not need to link with SFML, 
// it just needs some definitions for its helper functions.

#include <AntTweakBar.h>


//  TwEventSFML returns zero if msg has not been handled, 
//  and a non-zero value if it has been handled by the AntTweakBar library.
int TW_CALL TwEventSFML(const void *sfmlEvent, unsigned char majorVersion, unsigned char minorVersion)
{
    // Assume version 1.6 (will possibly not work for version != 1.6, but give it a chance)
    /*
    if (majorVersion > 1 || (majorVersion == 1 && minorVersion > 6)
    {
        static const char *g_ErrBadSFMLVersion = "Unsupported SFML version";
        TwSetLastError(g_ErrBadSFMLVersion);
        return 0;
    }
    */
    (void)majorVersion, (void)minorVersion;

    int handled = 0;
    const sf::Event *event = (const sf::Event *)sfmlEvent;
    TwMouseAction mouseAction;
    int key = 0;
    static int s_KMod = 0;
    static bool s_PreventTextHandling = false;
    static int s_WheelPos = 0;

    if (event == NULL)
        return 0;

    switch (event->Type)
    {
    case sf::Event::KeyPressed:
        s_PreventTextHandling = false;
        s_KMod = 0;
        if (event->Key.Shift)   s_KMod |= TW_KMOD_SHIFT;
        if (event->Key.Alt)     s_KMod |= TW_KMOD_ALT;
        if (event->Key.Control) s_KMod |= TW_KMOD_CTRL;
        key = 0;
        switch (event->Key.Code)
        {
        case sf::Key::Escape:
            key = TW_KEY_ESCAPE;
            break;
        case sf::Key::Return:
            key = TW_KEY_RETURN;
            break;
        case sf::Key::Tab:
            key = TW_KEY_TAB;
            break;
        case sf::Key::Back:
            key = TW_KEY_BACKSPACE;
            break;
        case sf::Key::PageUp:
            key = TW_KEY_PAGE_UP;
            break;
        case sf::Key::PageDown:
            key = TW_KEY_PAGE_DOWN;
            break;
        case sf::Key::Up:
            key = TW_KEY_UP;
            break;
        case sf::Key::Down:
            key = TW_KEY_DOWN;
            break;
        case sf::Key::Left:
            key = TW_KEY_LEFT;
            break;
        case sf::Key::Right:
            key = TW_KEY_RIGHT;
            break;
        case sf::Key::End:
            key = TW_KEY_END;
            break;
        case sf::Key::Home:
            key = TW_KEY_HOME;
            break;
        case sf::Key::Insert:
            key = TW_KEY_INSERT;
            break;
        case sf::Key::Delete:
            key = TW_KEY_DELETE;
            break;
        case sf::Key::Space:
            key = TW_KEY_SPACE;
            break;
        default:
            if (event->Key.Code >= sf::Key::F1 && event->Key.Code <= sf::Key::F15)
                key = TW_KEY_F1 + event->Key.Code - sf::Key::F1;
            else if (s_KMod & TW_KMOD_ALT) 
            {
                if (event->Key.Code >= sf::Key::A && event->Key.Code <= sf::Key::Z)
                {
                    if (s_KMod & TW_KMOD_SHIFT)
                        key = 'A' + event->Key.Code - sf::Key::A;
                    else
                        key = 'a' + event->Key.Code - sf::Key::A;
                }
            }
        }
        if (key != 0) 
        {
            handled = TwKeyPressed(key, s_KMod);
            s_PreventTextHandling = true;
        }
        break;
    case sf::Event::KeyReleased:
        s_PreventTextHandling = false;
        s_KMod = 0;
        break;
    case sf::Event::TextEntered:
        if (!s_PreventTextHandling && event->Text.Unicode != 0 && (event->Text.Unicode & 0xFF00) == 0)
        {
            if ((event->Text.Unicode & 0xFF) < 32) // CTRL+letter
                handled = TwKeyPressed((event->Text.Unicode & 0xFF)+'a'-1, TW_KMOD_CTRL|s_KMod);
            else 
                handled = TwKeyPressed(event->Text.Unicode & 0xFF, 0);
        }
        s_PreventTextHandling = false;
        break;
    case sf::Event::MouseMoved:
        handled = TwMouseMotion(event->MouseMove.X, event->MouseMove.Y);
        break;
    case sf::Event::MouseButtonPressed:
    case sf::Event::MouseButtonReleased:
        mouseAction = (event->Type==sf::Event::MouseButtonPressed) ? TW_MOUSE_PRESSED : TW_MOUSE_RELEASED;
        switch (event->MouseButton.Button) 
        {
        case sf::Mouse::Left:
            handled = TwMouseButton(mouseAction, TW_MOUSE_LEFT);
            break;
        case sf::Mouse::Middle:
            handled = TwMouseButton(mouseAction, TW_MOUSE_MIDDLE);
            break;
        case sf::Mouse::Right:
            handled = TwMouseButton(mouseAction, TW_MOUSE_RIGHT);
            break;
        default:
            break;
        }
        break;
    case sf::Event::MouseWheelMoved:
        s_WheelPos += event->MouseWheel.Delta;
        handled = TwMouseWheel(s_WheelPos);
        break;
    case sf::Event::Resized:
        // tell the new size to TweakBar
        TwWindowSize(event->Size.Width, event->Size.Height);
        // do not set 'handled', sf::Event::Resized may be also processed by the client application
        break;
    default:
        break;
    }

    return handled;
}
