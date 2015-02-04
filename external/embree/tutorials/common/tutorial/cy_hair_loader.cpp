#include "cy_hair_loader.h"

#include <stdio.h>
#include <math.h>

namespace embree
{
  // -------------------------------------------
  // -- Based on the format by Cem Yuksel     --
  // -- www.cemyuksel.com/research/hairmodels --
  // -------------------------------------------

#define CY_HAIR_FILE_SEGMENTS_BIT	1
#define CY_HAIR_FILE_POINTS_BIT		2
#define CY_HAIR_FILE_THICKNESS_BIT	4
#define CY_HAIR_FILE_TRANSPARENCY_BIT	8
#define CY_HAIR_FILE_COLORS_BIT		16
#define CY_HAIR_FILE_INFO_SIZE		88

  struct cyHeader
  {
    char		signature[4];	///< This should be "HAIR"
    unsigned int	numStrands;	///< number of hair strands
    unsigned int	numPoints;	///< total number of points of all strands
    unsigned int	bitarrays;		///< bit array of data in the file
    unsigned int	defaultSegments;      ///< default number of segments of each strand
    float		defaultThickness;	        ///< default thickness of hair strands
    float		defaultTransparency;	        ///< default transparency of hair strands
    float		defaultColor[3];		///< default color of hair strands
    char		info[CY_HAIR_FILE_INFO_SIZE];	///< information about the file
  };

  
  class cyHairFile
  {
  public:
    cyHairFile() : segments(NULL), points(NULL), thickness(NULL), transparency(NULL), colors(NULL) {init(); }
    ~cyHairFile() {
      if ( segments ) delete [] segments;
      if ( points ) delete [] points;
      if ( colors ) delete [] colors;
      if ( thickness ) delete [] thickness;
      if ( transparency ) delete [] transparency;      
    }

   
    void init()
    {
      header.signature[0] = 'H';
      header.signature[1] = 'A';
      header.signature[2] = 'I';
      header.signature[3] = 'R';
      header.numStrands = 0;
      header.numPoints = 0;
      header.bitarrays = 0;	
      header.defaultSegments = 0;
      header.defaultThickness = 1.0f;
      header.defaultTransparency = 0.0f;
      header.defaultColor[0] = 1.0f;
      header.defaultColor[1] = 1.0f;
      header.defaultColor[2] = 1.0f;
      memset( header.info, '\0', CY_HAIR_FILE_INFO_SIZE );
    }

    int load( const char *filename )
    {
      init();

      FILE *file;
      file = fopen( filename, "rb" );
      if ( file == NULL )
        THROW_RUNTIME_ERROR("can't open file");

      size_t h = fread( &header, sizeof(cyHeader), 1, file );

      if ( h < 1 ) 
        THROW_RUNTIME_ERROR("can't read header");

      if ( strncmp( header.signature, "HAIR", 4) != 0 ) 
        THROW_RUNTIME_ERROR("wrong signature");

      if ( header.bitarrays & CY_HAIR_FILE_SEGMENTS_BIT ) {
        segments = new unsigned short[ header.numStrands ];
        size_t r = fread( segments, sizeof(unsigned short), header.numStrands, file );
        if ( r < header.numStrands ) 
          THROW_RUNTIME_ERROR("error reading segments");
      }

      if ( header.bitarrays & CY_HAIR_FILE_POINTS_BIT ) {
        points = new float[ header.numPoints*3 ];
        size_t r = fread( points, sizeof(float), header.numPoints*3, file );
        if ( r < header.numPoints*3 ) 
          THROW_RUNTIME_ERROR("error reading points");
      }

      if ( header.bitarrays & CY_HAIR_FILE_THICKNESS_BIT ) {
        thickness = new float[ header.numPoints ];
        size_t r = fread( thickness, sizeof(float), header.numPoints, file );
        if ( r < header.numPoints ) 
          THROW_RUNTIME_ERROR("error reading thickness values");
      }

      if ( header.bitarrays & CY_HAIR_FILE_TRANSPARENCY_BIT ) {
        transparency = new float[ header.numPoints ];
        size_t r = fread( transparency, sizeof(float), header.numPoints, file );
        if ( r < header.numPoints ) 
          THROW_RUNTIME_ERROR("error reading transparency values");          
      }

      if ( header.bitarrays & CY_HAIR_FILE_COLORS_BIT ) {
        colors = new float[ header.numPoints*3 ];
        size_t r = fread( colors, sizeof(float), header.numPoints*3, file );
        if ( r < header.numPoints*3 ) 
          THROW_RUNTIME_ERROR("error reading color values");          
      }

      fclose( file );

      return header.numStrands;
    }


    cyHeader header;
    unsigned short	*segments;
    float		*points;
    float		*thickness;
    float		*transparency;
    float		*colors;

    
  };

  void loadCYHair(const FileName& fileName, OBJScene& scene, Vec3fa& offset)
  {
    cyHairFile cyFile;
    int numHairs = cyFile.load(fileName.c_str());
    std::cout << "Successfully loaded hair data" << std::endl; 

    OBJScene::HairSet* hairset = new OBJScene::HairSet; 

    for (size_t i=0;i<cyFile.header.numPoints;i++)
      {
        Vec3fa v;
        v.x = cyFile.points[3*i+0];
        v.y = cyFile.points[3*i+1];
        v.z = cyFile.points[3*i+2];
        v.w = cyFile.thickness ? cyFile.thickness[i] : 0.1f;
        v.x-=offset.x;
        v.y-=offset.y;
        v.z-=offset.z;
        hairset->v.push_back(v);
      }

    ssize_t index = 0;
    for (ssize_t i=0;i<cyFile.header.numStrands;i++)
      {
        if (cyFile.segments)
          {
            ssize_t numSegments = cyFile.segments[i];       
            for (ssize_t j=0; j<numSegments-3; j+=3)
              {
                hairset->hairs.push_back(OBJScene::Hair(index + j,i));
              }
            index += numSegments+1;	
          }
        else
          {
            ssize_t numSegments = cyFile.header.defaultSegments;       
            for (ssize_t j=0; j<numSegments-3; j+=3)
              {
                hairset->hairs.push_back(OBJScene::Hair(index + j,i));
              }
            index += numSegments+1;	

          }
      }

    scene.hairsets.push_back(hairset);
  }
}
