#include "YImage.hpp"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstring>

#include <png.h>
#include <zlib.h> // For Z_BEST_SPEED

// for: #define offsetof(TYPE, MEMBER) ((size_t) &((TYPE *)0)->MEMBER)
// we use this to determine the pixel format on the fly
#include <cstddef>

YImage::YImage()
	: m_width( 0 ), m_height( 0 ), m_data( NULL )
{
	assert( sizeof( YPixel ) == 4 && "YPixel struct shouldn't be padded" ) ;
}

YImage::YImage( const YImage& rhs )
	: m_width( 0 ), m_height( 0 ), m_data( NULL )
{
	*this = rhs ;
}

YImage&
YImage::operator=( const YImage& rhs )
{
    if( !rhs.m_data )
    {
        assert( 0 == rhs.m_width ) ;
        assert( 0 == rhs.m_height ) ;
        
        if( m_data ) free( m_data ) ;
        m_data = NULL ;
        m_width = 0 ;
        m_height = 0 ;
        
        return *this ;
    }
    
	if( m_width != rhs.m_width || m_height != rhs.m_height )
	{
		if( m_data ) free( m_data ) ;
		
		m_width = rhs.m_width ;
		m_height = rhs.m_height ;
		m_data = (YPixel*) malloc( m_width * m_height * sizeof(YPixel) ) ;
	}
	
	assert( m_data ) ;
	memcpy( m_data, rhs.m_data, m_width * m_height * sizeof(YPixel) ) ;
	
	return *this ;
}

YImage::~YImage()
{
	if( m_data ) free( m_data ) ;
	m_data = NULL ;
	m_width = m_height = 0 ;
}

YImage::YPixel*
YImage::data()
{
	return m_data ;
}
const YImage::YPixel*
YImage::data()
const
{
	return m_data ;
}
YImage::YPixel&
YImage::at( int i, int j )
{
	assert( i >= 0 && i < m_width ) ;
	assert( j >= 0 && j < m_height ) ;
	
	return m_data[ i + j * m_width ] ;
}
const YImage::YPixel&
YImage::at( int i, int j )
const
{
	assert( i >= 0 && i < m_width ) ;
	assert( j >= 0 && j < m_height ) ;
	
	return m_data[ i + j * m_width ] ;
}

int YImage::width() const
{
	return m_width ;
}
int YImage::height() const
{
	return m_height ;
}

// Creates a new image with the specified size, preserving as
// much as possible the old image (as in a window resize).
void YImage::resize( int widthh, int heightt )
{
	// do nothing
	if( m_width == widthh && m_height == heightt )
	{
		assert( m_data ) ;
		return ;
	}
	
	YPixel* new_data = (YPixel*) malloc( widthh * heightt * sizeof(YPixel) ) ;
	memset( new_data, 0, widthh * heightt * sizeof(YPixel) ) ;
	
	if( m_data )
	{
		int min_width = std::min( m_width, widthh ) ;
		int min_height = std::min( m_height, heightt ) ;
		
		for( int j = 0 ; j < min_height ; ++j )
			memcpy( new_data + j*widthh, m_data + j*m_width, min_width * sizeof(YPixel) ) ;
		
		free( m_data ) ;
	}
	
	m_width = widthh ;
	m_height = heightt ;
	m_data = new_data ;
}


void YImage::greyscale()
{
	for( int i = 0 ; i < m_width * m_height ; ++i )
	{
		YPixel* pix = m_data + i ;
		int greyval = (int) pix->r + (int) pix->g + (int) pix->b ;
		greyval /= 3 ;
		pix->r = pix->g = pix->b = greyval ;
	}
}

void YImage::flip()
{
	for( int j = 0 ; j < m_height / 2 ; ++j )
	for( int i = 0 ; i < m_width ; ++i )
	{
		YPixel* lhs = &at( i,j ) ;
		YPixel* rhs = &at( i, m_height - j - 1 ) ;
		
		YPixel temp = *lhs ;
		
		*lhs = *rhs ;
		*rhs = temp ;
	}
}

void YImage::mirror()
{
	for( int j = 0 ; j < m_height ; ++j )
	for( int i = 0 ; i < m_width / 2 ; ++i )
	{
		YPixel* lhs = &at( i,j ) ;
		YPixel* rhs = &at( m_width - i - 1, j ) ;
		
		YPixel temp = *lhs ;
		
		*lhs = *rhs ;
		*rhs = temp ;
	}
}

bool YImage::same( const YImage& rhs ) const
{
    if( m_height != rhs.m_height || m_width != rhs.m_width ) return false;
    
    // Check if this->m_data and rhs.m_data point to the same data.
    // This also checks for if both are NULL.
    if( m_data == rhs.m_data ) return true;
    
    assert( m_data && rhs.m_data );
    return 0 == memcmp( m_data, rhs.m_data, m_width * m_height * sizeof(YPixel) );
}
bool YImage::same_rgb( const YImage& rhs ) const
{
    if( m_height != rhs.m_height || m_width != rhs.m_width ) return false;
    
    // Check if this->m_data and rhs.m_data point to the same data.
    // This also checks for if both are NULL.
    if( m_data == rhs.m_data ) return true;
    
    assert( m_data && rhs.m_data );
    
    for( int i = 0; i < m_height*m_width; ++i )
    {
        if(
            m_data[i].r != rhs.m_data[i].r ||
            m_data[i].g != rhs.m_data[i].g ||
            m_data[i].b != rhs.m_data[i].b )
            return false;
    }
    
    return true;
}


// We use out own reading/writing functions because libpng may have been compiled using a
// different compiler & libc.  That could make the FILE*'s incompatible.
static void user_write_data(png_structp png_ptr, png_bytep data, png_size_t length);
static void user_write_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
	png_voidp write_io_ptr = png_get_io_ptr(png_ptr) ;
	fwrite( (unsigned char*) data, length, 1, (FILE*) write_io_ptr ) ;
}

static void user_flush_data(png_structp png_ptr);
static void user_flush_data(png_structp png_ptr)
{
	png_voidp write_io_ptr = png_get_io_ptr(png_ptr) ;
	fflush( (FILE*) write_io_ptr ) ;
}

static void user_read_data(png_structp png_ptr, png_bytep data, png_size_t length) ;
static void user_read_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
	png_voidp read_io_ptr = png_get_io_ptr(png_ptr) ;
	(void) fread( (unsigned char*) data, length, 1, (FILE*) read_io_ptr ) ;
}

// Some error strings
#define ERROR_STRING_OPEN "Error opening %s\n"
#define ERROR_STRING_WRITING "libpng encountered a problem writing %s\n"
#define ERROR_STRING_INVALID_FILE "Invalid png file: %s\n"
#define ERROR_STRING_LIBERR "libpng encountered a problem (not related to file)\n"
#define ERROR_STRING_MAYBE_FILE "libpng encountered a problem (may be related to file: %s)\n"

bool YImage::save( const char* fname, const bool fast )
const
{
	FILE* fp = NULL ;
	bool rval = true ;
	png_structp png_ptr = NULL ;
	png_infop info_ptr = NULL ;
	
	// Open the file for reading in binary mode.
	fp = fopen( fname, "wb" ) ;
	if( !fp )
	{
		fprintf( stderr, ERROR_STRING_OPEN, fname ) ;
		rval = false ;
		goto YImage_save_cleanup ;
	}
	
	// Allocate the png structs.
	png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL ) ;
	if( !png_ptr )
	{
		fprintf( stderr, ERROR_STRING_WRITING, fname ) ;
		rval = false ;
		goto YImage_save_cleanup ;
	}
	info_ptr = png_create_info_struct( png_ptr ) ;
	if( !info_ptr )
	{
		fprintf( stderr, ERROR_STRING_WRITING, fname ) ;
		rval = false ;
		goto YImage_save_cleanup ;
    }
	
	// Set up the png error routine.
	if( setjmp(png_jmpbuf(png_ptr)) )
	{
		fprintf( stderr, ERROR_STRING_WRITING, fname ) ;
		rval = false ;
		goto YImage_save_cleanup ;
    }
	
	// Give libpng the FILE*.
	// png_init_io( png_ptr, fp ) ;
	// or
	// use our own write callback
	png_set_write_fn( png_ptr, fp, (png_rw_ptr) user_write_data, user_flush_data ) ;
	
	// We'll use the low-level interface since the high-level interface won't handle
	// png_set_filler() which we need to tell libpng to strip out the A from our
	// 4-byte pixels.
	
	// First we set and write the info struct.
	png_set_IHDR(
		png_ptr, info_ptr,
		m_width, m_height,
		8, PNG_COLOR_TYPE_RGB_ALPHA,
		PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT
		) ;
	png_write_info( png_ptr, info_ptr ) ;
	
	// If we've been asked to write quickly, speed up the compression.
	if( fast ) png_set_compression_level( png_ptr, Z_BEST_SPEED );
	
	// Then we set up transforms.
	/*
	// 1. tell libpng to strip out the filler byte in our 4-byte pixels
	// if YPixel::a comes after any other member (b,g,r), we strip AFTER
	if( offsetof( YPixel, a ) > offsetof( YPixel, b ) ) {
		png_set_filler( png_ptr, 0, PNG_FILLER_AFTER ) ;
		// printf("alpha after\n");
	} else {
		png_set_filler( png_ptr, 0, PNG_FILLER_BEFORE ) ;
		// printf("alpha before\n");
	}
	*/
	if( offsetof( YPixel, a ) < offsetof( YPixel, b ) )
		png_set_swap_alpha( png_ptr ) ;
	// 2. tell libpng how our color triples are stored (b < r or vice versa)
	if( offsetof( YPixel, b ) < offsetof( YPixel, r ) )
	{
		png_set_bgr( png_ptr ) ;
		// printf("bgr\n") ;
	}
	// else { printf("rgb\n"); }
	// printf( "offsetof r, b: %d %d\n", offsetof( YPixel, r ), offsetof( YPixel, b ) );
	
	// Finally we create a row_pointers[] pointing into our data* and write the png out to the FILE*.
	{
		// 1. allocate row pointers array
		png_bytep* row_pointers = (png_bytep*) png_malloc( png_ptr, m_height * sizeof(png_bytep) ) ;
		// 2. point row pointers into m_data
		for( int i = 0 ; i < m_height ; ++i ) {
			row_pointers[i] = (png_bytep) (m_data + i*m_width) ;
		}
		// 3. write the image data
		png_write_image( png_ptr, row_pointers ) ;
		// 4. free row pointers array
		png_free( png_ptr, row_pointers ) ;
	}
	
	// Write out end info.  We're done.  Fall through to cleanup.
	png_write_end( png_ptr, NULL ) ;
	
YImage_save_cleanup:
	png_destroy_write_struct( png_ptr ? &png_ptr : NULL, info_ptr ? &info_ptr : NULL ) ;
	if( fp ) fclose( fp ) ;
	
	return rval ;
}

bool YImage::load(const char* fname)
{
	FILE* fp = NULL ;
	bool rval = true;
	png_structp png_ptr = NULL ;
	png_infop info_ptr = NULL ;
    
	// for checking the png header
	const size_t PNG_BYTES_TO_CHECK = 4 ; // example.c uses 4
	png_byte header[ PNG_BYTES_TO_CHECK ] ;
	
	// Open the file for reading in binary mode.
	fp = fopen( fname, "rb" ) ;
	if( !fp )
	{
		fprintf( stderr, ERROR_STRING_OPEN, fname ) ;
		rval = false;
		goto YImage_load_cleanup ;
	}
	
	// Check some bytes at the beginning of the file to make sure it's a png.
	if( PNG_BYTES_TO_CHECK != fread( header, 1, PNG_BYTES_TO_CHECK, fp ) )
	{
		fprintf( stderr, ERROR_STRING_INVALID_FILE, fname ) ;
		rval = false;
		goto YImage_load_cleanup ;
	}
	if( png_sig_cmp( header, 0, PNG_BYTES_TO_CHECK ) )
	{
		fprintf( stderr, ERROR_STRING_INVALID_FILE, fname ) ;
		rval = false;
		goto YImage_load_cleanup ;
	}
	
	// Since it looks like we have a good png file, allocate the png structs.
	png_ptr = png_create_read_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL ) ;
	if( !png_ptr )
	{
		fprintf( stderr, ERROR_STRING_LIBERR ) ;
		rval = false;
		goto YImage_load_cleanup ;
	}
	info_ptr = png_create_info_struct( png_ptr ) ;
	if( !info_ptr )
	{
		fprintf( stderr, ERROR_STRING_LIBERR ) ;
		rval = false;
		goto YImage_load_cleanup ;
    }
	
	// Set up the png error routine.
	if( setjmp(png_jmpbuf(png_ptr)) )
	{
		fprintf( stderr, ERROR_STRING_MAYBE_FILE, fname ) ;
		rval = false;
		goto YImage_load_cleanup ;
    }
	
	// Give libpng the FILE*, tell it how many bytes we read.
	// png_init_io( png_ptr, fp ) ;
	// or
	// use our own read callback
	png_set_read_fn( png_ptr, fp, (png_rw_ptr) user_read_data ) ;
	
	png_set_sig_bytes( png_ptr, PNG_BYTES_TO_CHECK ) ;
	
	// We'll use the low-level interface since the high-level interface won't handle
	// png_set_filler() which we need to guarantee there'll be a filler "A" in our
	// 4-byte ARGB pixels.
	// Really the low-level interface isn't more complicated than the high-level interface.
	// To choose transform flags we have to query the png_info struct.
	// Instead of OR'ing in another transform flag (high-level interface), we call a set
	// method (low-level interface).
	
	// First we read the info struct.
	png_read_info( png_ptr, info_ptr ) ;
	
	// Now we set up transforms.
	// 1. convert gray and paletted to rgb (this guarantees 8 or 16 bit depths (rgb must be 8 or 16)).
	//    also expand the alpha color to an alpha channel and convert 16 bit depths to 8.
	// 2. if we don't have alpha, add an opaque channel.  also swap RGB and BGR depending on YPixel.
	// 3. ask libpng to deinterlace
	{
		
		png_byte color_type = png_get_color_type( png_ptr, info_ptr ) ;
		png_byte depth = png_get_bit_depth( png_ptr, info_ptr ) ;
		
		// 1
		if( color_type == PNG_COLOR_TYPE_PALETTE )
			png_set_expand( png_ptr ) ;
		if( color_type == PNG_COLOR_TYPE_GRAY && depth < 8 )
			png_set_expand( png_ptr ) ;
		if( png_get_valid( png_ptr, info_ptr, PNG_INFO_tRNS ) )
			png_set_expand( png_ptr ) ;
		if( depth == 16 )
			png_set_strip_16( png_ptr ) ;
		if( color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA )
			png_set_gray_to_rgb( png_ptr ) ;
		
		
		// NOTE: This next step affects the layout of the channels in the pixel.
		//       We turn the pixel into the format in YPixel,
		//       with the restrication that alpha comes at the beginning or end
		//       and green is sandwiched between red and blue.
		//       The possibilities are: ARGB, ABGR, RGBA, BGRA.
		
		// 2
		if( color_type != PNG_COLOR_TYPE_GRAY_ALPHA && color_type != PNG_COLOR_TYPE_RGB_ALPHA && !png_get_valid( png_ptr, info_ptr, PNG_INFO_tRNS ) )
			png_set_filler( png_ptr, 0xFF, offsetof( YPixel, a ) < offsetof( YPixel, b ) ? PNG_FILLER_BEFORE : PNG_FILLER_AFTER ) ;
		else if( offsetof( YPixel, a ) < offsetof( YPixel, b ) )
				png_set_swap_alpha( png_ptr ) ;
		
		if( offsetof( YPixel, b ) < offsetof( YPixel, r ) )
			png_set_bgr( png_ptr ) ;
		
		// 3
		png_set_interlace_handling( png_ptr ) ;
		
	}
	
	// We're almost ready to copy over the png data.
	// First we must resize our data* and set our width & height.
	{
		png_uint_32 widthh = png_get_image_width( png_ptr, info_ptr ) ;
		png_uint_32 heightt = png_get_image_height( png_ptr, info_ptr ) ;
		
		// fprintf( stderr, "width: %d, height: %d\n", (int) widthh, (int) heightt ) ;
		
		resize( widthh, heightt ) ;
	}
	
	// Now we can create a rows[] pointing into our data* and read the png into our buffer.
	{
		// fprintf( stderr, "width: %d, height: %d\n", (int) m_width, (int) m_height ) ;
		/*
		
		png_byte channels = png_get_channels( png_ptr, info_ptr ) ;
		assert( 4 == channels ) ;
		png_byte depth = png_get_bit_depth( png_ptr, info_ptr ) ;
		assert( 8 == depth ) ;
		png_uint_32 rowbytes = png_get_rowbytes( png_ptr, info_ptr ) ;
		assert( sizeof(YPixel) * m_width == rowbytes ) ;
		*/
		
		// 1. allocate row pointers array
		png_bytep* row_pointers = (png_bytep*) png_malloc( png_ptr, m_height * sizeof(png_bytep) ) ;
		// 2. point row pointers into m_data
		for( int i = 0 ; i < m_height ; ++i )
			row_pointers[i] = (png_bytep) (m_data + i*m_width ) ;
		// 3. read the image data
		png_read_image( png_ptr, row_pointers ) ;
		// 4. free row pointers array
		png_free( png_ptr, row_pointers ) ;
	}
	
	// Read the end info.  We're done.  Fall through to cleanup.
	png_read_end( png_ptr, NULL ) ;
	
YImage_load_cleanup:
	// due to (what looks like) a bug in libpng-1.2.4, we can't pass NULL for the png_ptr arg
	if( png_ptr ) png_destroy_read_struct( &png_ptr, info_ptr ? &info_ptr : NULL, NULL ) ;
	
	if( fp ) fclose( fp ) ;
	
	return rval ;
}
