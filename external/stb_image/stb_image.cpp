// dummy cpp file to compile the code in case IGL is in header only mode
#ifdef STBI_LOAD
	#define STB_IMAGE_IMPLEMENTATION
	#include "stb_image.h"
#endif

#ifdef STBI_WRITE
	#define STB_IMAGE_WRITE_IMPLEMENTATION
	#include "stb_image_write.h"
#endif

 