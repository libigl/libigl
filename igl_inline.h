// This should *NOT* be contained in a IGL_*_H ifdef, since it may be defined
// differently based on when it is included
#ifdef IGL_INLINE
#undef IGL_INLINE
#endif

#ifdef IGL_HEADER_ONLY
#  define IGL_INLINE inline
#else
#  define IGL_INLINE
#endif
