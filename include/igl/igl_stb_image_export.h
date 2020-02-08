
#ifndef IGL_STB_IMAGE_EXPORT_H
#define IGL_STB_IMAGE_EXPORT_H

#ifdef IGL_STB_IMAGE_STATIC_DEFINE
#  define IGL_STB_IMAGE_EXPORT
#  define IGL_STB_IMAGE_NO_EXPORT
#else
#  ifndef IGL_STB_IMAGE_EXPORT
#    ifdef igl_stb_image_EXPORTS
        /* We are building this library */
#      define IGL_STB_IMAGE_EXPORT 
#    else
        /* We are using this library */
#      define IGL_STB_IMAGE_EXPORT 
#    endif
#  endif

#  ifndef IGL_STB_IMAGE_NO_EXPORT
#    define IGL_STB_IMAGE_NO_EXPORT 
#  endif
#endif

#ifndef IGL_STB_IMAGE_DEPRECATED
#  define IGL_STB_IMAGE_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef IGL_STB_IMAGE_DEPRECATED_EXPORT
#  define IGL_STB_IMAGE_DEPRECATED_EXPORT IGL_STB_IMAGE_EXPORT IGL_STB_IMAGE_DEPRECATED
#endif

#ifndef IGL_STB_IMAGE_DEPRECATED_NO_EXPORT
#  define IGL_STB_IMAGE_DEPRECATED_NO_EXPORT IGL_STB_IMAGE_NO_EXPORT IGL_STB_IMAGE_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef IGL_STB_IMAGE_NO_DEPRECATED
#    define IGL_STB_IMAGE_NO_DEPRECATED
#  endif
#endif

#endif /* IGL_STB_IMAGE_EXPORT_H */
