// Adapted from QuickLook-pfm-master/src/GenerateThumbnailForURL.m by Andreas
// Steinel (https://github.com/lnxbil/quicklook-pfm)
#import <CoreFoundation/CoreFoundation.h>
#import <CoreServices/CoreServices.h>
#import <QuickLook/QuickLook.h>
#import <Foundation/Foundation.h>
#import <AppKit/NSImage.h>

#include "render_to_buffer.h"

// Creates a bitmap context for ARGB images (not sure if we really need a
// separate function for this... I think there should be something builtin that
// works just as well)
// Taken from https://github.com/lnxbil/quicklook-pfm
CGContextRef CreateARGBBitmapContext (CGSize size)
{
    CGContextRef    context = NULL;
    CGColorSpaceRef colorSpace;
    void *          bitmapData;
    int             bitmapByteCount;
    int             bitmapBytesPerRow;
    
    // Get image width, height. We'll use the entire image.
    size_t pixelsWide = (size_t) size.width; //CGImageGetWidth(inImage);
    size_t pixelsHigh = (size_t) size.height; //CGImageGetHeight(inImage);
    
    // Declare the number of bytes per row. Each pixel in the bitmap in this
    // example is represented by 4 bytes; 8 bits each of red, green, blue, and
    // alpha.
    bitmapBytesPerRow   = (pixelsWide * 4);
    bitmapByteCount     = (bitmapBytesPerRow * pixelsHigh);
    
    // Use the generic RGB color space.
    // FIXME: Das verstehe ich net. Die Doku sagt nicht, dass es Deprecated ist und
    //        gibt auch keine Auswahlmöglichkeit an :-(
    colorSpace = CGColorSpaceCreateWithName(kCGColorSpaceGenericRGB);
    if (colorSpace == NULL)
    {
        NSLog(@"Error allocating color space\n");
        return NULL;
    }
    
    // Allocate memory for image data. This is the destination in memory
    // where any drawing to the bitmap context will be rendered.
    bitmapData = malloc( bitmapByteCount );
    if (bitmapData == NULL) 
    {
        NSLog(@"Memory not allocated!");
        CGColorSpaceRelease( colorSpace );
        return NULL;
    }
    
    // Create the bitmap context. We want pre-multiplied ARGB, 8-bits 
    // per component. Regardless of what the source image format is 
    // (CMYK, Grayscale, and so on) it will be converted over to the format
    // specified here by CGBitmapContextCreate.
    context = CGBitmapContextCreate (bitmapData,
                                     pixelsWide,
                                     pixelsHigh,
                                     8,      // bits per component
                                     bitmapBytesPerRow,
                                     colorSpace,
                                     kCGBitmapByteOrderDefault|kCGImageAlphaPremultipliedLast);
    if (context == NULL)
    {
        free (bitmapData);
        NSLog(@"Context not created!");
    }
    
    // Make sure and release colorspace before returning
    CGColorSpaceRelease( colorSpace );
    
    return context;
}




// Generate a thumbnail for a given file
OSStatus GenerateThumbnailForURL(
  void *thisInterface, 
  QLThumbnailRequestRef thumbnail, 
  CFURLRef url, 
  CFStringRef contentTypeUTI, 
  CFDictionaryRef options, 
  CGSize maxSize)
{
    if (QLThumbnailRequestIsCancelled(thumbnail))
        return noErr;
    
    NSAutoreleasePool* pool = [[NSAutoreleasePool alloc] init];
    
    CFStringRef file = CFURLCopyFileSystemPath(url,kCFURLPOSIXPathStyle);
    
    if (QLThumbnailRequestIsCancelled(thumbnail))
        return noErr;

    // raw pixel data memory of 64 * 64 pixel size
    // http://stackoverflow.com/a/3797923/148668
    // Multi-sampling for anti-aliasing
    const double MS = 2;
    const int width = maxSize.width;
    // Widescreen aspect ratio
    const int height = (270./460.)*maxSize.height;
    const int bwidth = MS*width;
    const int bheight = MS*height;
    const int bchannels = 4;
    UInt8 buffer[bwidth * bheight * bchannels];
    const char *cs = CFStringGetCStringPtr( file, CFStringGetSystemEncoding()) ;
    char cs_buf[1024];
    // http://stackoverflow.com/a/1609664/148668
    if(cs == NULL)
    {
      CFStringGetCString(file, cs_buf, [(NSString*)file length]+1, kCFStringEncodingUTF8);
      cs = cs_buf;
    }
    const float OPAQUE_WHITE[4] = {1,1,1,1};
    render_to_buffer(cs,OPAQUE_WHITE,bwidth,bheight,buffer);
    CGColorSpaceRef colorspace = CGColorSpaceCreateDeviceRGB();
    CFDataRef rgbData = CFDataCreate(NULL, buffer, bwidth * bheight * 4);
    CGDataProviderRef provider = CGDataProviderCreateWithCFData(rgbData);
    // http://forum.sparrow-framework.org/topic/create-uiimage-from-pixel-data-problems
    CGImageRef image = CGImageCreate(
      bwidth,
      bheight,
      8,
      8*4,
      bwidth * 4,
      colorspace,
      kCGBitmapByteOrderDefault|kCGImageAlphaPremultipliedLast,
      provider,
      NULL,
      true,
      kCGRenderingIntentDefault); 
    CFRelease(rgbData);
    CGDataProviderRelease(provider);
    CGColorSpaceRelease(colorspace);
    NSSize size = {width,height};
    
    if (QLThumbnailRequestIsCancelled(thumbnail))
        return noErr;

    // Draw onto context as textured rectangle
    CGContextRef cgctx = CreateARGBBitmapContext(size);
    CGRect rect = CGRectMake(0,0, width, height);
    CGContextClearRect(cgctx,rect);
    CGContextDrawImage(cgctx, rect, image); 
    CGImageRef newCGImage = CGBitmapContextCreateImage(cgctx);
    CGContextRelease(cgctx);
    
    QLThumbnailRequestSetImage(thumbnail, newCGImage, NULL);
    
    // Releasing image
    CGImageRelease(newCGImage);
    [pool release];
    return noErr;
}


void CancelThumbnailGeneration(void* thisInterface, QLThumbnailRequestRef thumbnail)
{
    // implement only if supported
}

