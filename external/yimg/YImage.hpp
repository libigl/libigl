#ifndef __YImage_hpp__
#define __YImage_hpp__

class YImage
{
public:
    // file loading/saving automatically picks up changes to this struct.
    // The possibilities are: ARGB, ABGR, RGBA, BGRA.
    struct YPixel {
        unsigned char r ;
        unsigned char g ;
        unsigned char b ;
        unsigned char a ;
    };
    
    YImage() ;
	YImage( const YImage& ) ;
	virtual ~YImage() ;
	
	YImage& operator=( const YImage& ) ;
	
	// Returns true if the images are the same size and have the same data (for r,g,b,a).
	bool same( const YImage& rhs ) const ;
	// Returns true if the images are the same size and have the same data for r,g,b.
	bool same_rgb( const YImage& rhs ) const ;
	
	// Saves a PNG image.  If 'fast' is true, sacrifices compression ratio for speed.
	bool save( const char* fname, const bool fast = false ) const ;
	bool load( const char* fname ) ;
	
	YPixel* data() ;
	const YPixel* data() const ;
	
	YPixel& at( int i, int j ) ;
	const YPixel& at( int i, int j ) const ;
	
	int width() const ;
	int height() const ;
	void resize( int width, int height ) ;
	
	// flip vertically
	void flip() ;
	// flip horizontally
	void mirror() ;
	// average rgb
	void greyscale() ;
	
protected:
	int m_width ;
	int m_height ;
	YPixel* m_data ; // raw image data
};

#endif /* __YImage_hpp__ */
