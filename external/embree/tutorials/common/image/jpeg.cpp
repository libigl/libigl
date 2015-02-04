// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#ifdef USE_LIBJPEG

#include "image/image.h"
#include "jpeglib.h"

namespace embree
{

    void compress(struct jpeg_compress_struct *cinfo, unsigned char *image)
    {

        /*! Start compression. */
        jpeg_start_compress(cinfo, TRUE);

        /*! Pointer to and size of a scanline in the image. */
        JSAMPROW scanline[1];  size_t bytes = cinfo->image_width * cinfo->input_components;

        /*! Here we use the library state variable 'next_scanline' as the loop index. */
        while (cinfo->next_scanline < cinfo->image_height) scanline[0] = &image[cinfo->next_scanline * bytes], jpeg_write_scanlines(cinfo, scanline, 1);

        /*! Finish compression. */
        jpeg_finish_compress(cinfo);

    }

    unsigned char *decompress(struct jpeg_decompress_struct *cinfo)
    {

        /*! Start decompression. */
        jpeg_start_decompress(cinfo);

        /*! Bytes per row in the scanline buffer. */
        size_t bytes = cinfo->output_width * cinfo->output_components;

        /*! Allocate scratch space for a single scanline. */
        JSAMPARRAY scanline = (*cinfo->mem->alloc_sarray)((j_common_ptr) cinfo, JPOOL_IMAGE, bytes, 1);

        /*! Allocate storage for the decompressed image. */
        unsigned char *image = (unsigned char *) malloc(cinfo->output_height * bytes);  if (image == NULL) return(NULL);

        /*! Here we use the library state variable 'output_scanline' as the loop index. */
        while (cinfo->output_scanline < cinfo->output_height) jpeg_read_scanlines(cinfo, scanline, 1), memcpy(&image[(cinfo->output_scanline - 1) * bytes], scanline[0], bytes);

        /*! Finish decompression. */
        jpeg_finish_decompress(cinfo);  return(image);

    }

    void encodeRGB8_to_JPEG(unsigned char *image, size_t width, size_t height, unsigned char **encoded, size_t *capacity)
    {

#if JPEG_LIB_VERSION >= 80

        /*! Compression parameters and scratch space pointers (allocated by the library). */
        struct jpeg_compress_struct cinfo;

        /*! The library error handler. */
        struct jpeg_error_mgr jerror;  cinfo.err = jpeg_std_error(&jerror);

        /*! Initialize the JPEG compression object. */
        jpeg_create_compress(&cinfo);

        /*! Specify the incoming image resolution, color space, and color space components. */
        cinfo.image_width = width;  cinfo.image_height = height;  cinfo.in_color_space = JCS_RGB;  cinfo.input_components = 3;

        /*! Fill in a sensible set of defaults. */
        jpeg_set_defaults(&cinfo);

        /*! Set the image quality. */
        jpeg_set_quality(&cinfo, 90, TRUE);

        /*! Specify the data source. */
        jpeg_mem_dest(&cinfo, encoded, capacity);

        /*! Compress and write the image into the target buffer. */
        compress(&cinfo, image);

        /*! At this point 'jerror.num_warnings' could be checked for corrupt-data warnings. */
        jpeg_destroy_compress(&cinfo);

#else  // JPEG_LIB_VERSION

        THROW_RUNTIME_ERROR("JPEG encoding into a memory buffer requires LibJPEG 8a or higher");

#endif // JPEG_LIB_VERSION

    }

    Ref<Image> loadJPEG(const FileName &filename)
    {

        /*! Open the source JPEG file. */
        FILE *file = fopen(filename.c_str(), "rb");  if (!file) THROW_RUNTIME_ERROR("Unable to open \"" + filename.str() + "\".");

        /*! Decompression parameters and scratch space pointers allocated by the library. */
        struct jpeg_decompress_struct cinfo;

        /*! The library error handler. */
        struct jpeg_error_mgr jerror;  cinfo.err = jpeg_std_error(&jerror);

        /*! Initialize the JPEG decompression object. */
        jpeg_create_decompress(&cinfo);

        /*! Specify the data source. */
        jpeg_stdio_src(&cinfo, file);

        /*! Read file parameters with jpeg_read_header(). */
        jpeg_read_header(&cinfo, TRUE);

        /*! Specify the color space and color space components of the decompressed image. */
        cinfo.out_color_space = JCS_RGB;  cinfo.output_components = 3;

        /*! Decompress the image into an output buffer and get the image dimensions. */
        unsigned char *rgb = decompress(&cinfo);  size_t width = cinfo.output_width;  size_t height = cinfo.output_height;

        /*! Allocate the Embree image. */
        Ref<Image> image = new Image4c(width, height, filename);

        /*! Convert the image from unsigned char RGB to unsigned char RGBA. */
        for (size_t y=0, i=0 ; y < height ; y++) {
          for (size_t x=0 ; x < width ; x++) {
            const float r = (float) rgb[i++] / 255.0f;
            const float g = (float) rgb[i++] / 255.0f;
            const float b = (float) rgb[i++] / 255.0f;
            image->set(x, y, Color4(r,g,b,1.0f));
          }
        }

        /*! Clean up. */
        jpeg_destroy_decompress(&cinfo);  free(rgb);  fclose(file);  return(image);

    }

    void storeJPEG(const Ref<Image> &image, const FileName &filename)
    {

        /*! Open the target JPEG file. */
        FILE *file = fopen(filename.c_str(), "wb");  if (!file) THROW_RUNTIME_ERROR("Unable to open \"" + filename.str() + "\".");

        /*! Compression parameters and scratch space pointers (allocated by the library). */
        struct jpeg_compress_struct cinfo;

        /*! The library error handler. */
        struct jpeg_error_mgr jerror;  cinfo.err = jpeg_std_error(&jerror);

        /*! Initialize the JPEG compression object. */
        jpeg_create_compress(&cinfo);

        /*! Specify the incoming image resolution, color space, and color space components. */
        cinfo.image_width = image->width;  cinfo.image_height = image->height;  cinfo.in_color_space = JCS_RGB;  cinfo.input_components = 3;

        /*! Fill in a sensible set of defaults. */
        jpeg_set_defaults(&cinfo);

        /*! Specify the data source. */
        jpeg_stdio_dest(&cinfo, file);

        /*! Allocate storage for the uncompressed packed image. */
        unsigned char *rgb = (unsigned char *) malloc(image->height * image->width);

        /*! Convert the image to unsigned char RGB. */
        for (size_t y=0, i=0 ; y < image->height ; y++) {
            for (size_t x=0 ; x < image->width ; x++) {
                const Color4 pixel = image->get(x, y);
                rgb[i++] = (unsigned char)(clamp(pixel.r) * 255.0f);
                rgb[i++] = (unsigned char)(clamp(pixel.g) * 255.0f);
                rgb[i++] = (unsigned char)(clamp(pixel.b) * 255.0f);
            }
        }

        /*! Compress and write the image into the target file. */
        compress(&cinfo, rgb);

        /*! At this point 'jerror.num_warnings' could be checked for corrupt-data warnings. */
        jpeg_destroy_compress(&cinfo);  free(rgb);  fclose(file);

    }

}

#endif // USE_LIBJPEG

