// Adapted from QuickLook-pfm-master/src/GeneratePreviewForURL.m by Andreas
// Steinel (https://github.com/lnxbil/quicklook-pfm)
#import <CoreFoundation/CoreFoundation.h>
#import <CoreServices/CoreServices.h>
#import <QuickLook/QuickLook.h>
#import <Foundation/Foundation.h>
#import <AppKit/NSImage.h>

#include "render_to_buffer.h"

// Generate a preview of a given file.
OSStatus GeneratePreviewForURL(
  void *thisInterface, 
  QLPreviewRequestRef preview, 
  CFURLRef url, 
  CFStringRef contentTypeUTI, 
  CFDictionaryRef options)
{
    if (QLPreviewRequestIsCancelled(preview))
        return noErr;
    
    NSAutoreleasePool* pool = [[NSAutoreleasePool alloc] init];
    
    CFStringRef file = CFURLCopyFileSystemPath(url,kCFURLPOSIXPathStyle);
    if (QLPreviewRequestIsCancelled(preview))
        return noErr;
    // http://stackoverflow.com/a/3797923/148668
    // Multi-sampling for anti-aliasing.
    const double MS = 2;
    const int width = 720;
    const int height = 405;
    const int bwidth = MS*width;
    const int bheight = MS*height;
    const int bchannels = 4;
    UInt8 buffer[bwidth * bheight * bchannels];

    const char *cs = CFStringGetCStringPtr( file, CFStringGetSystemEncoding()) ;
    char cs_buf[1024];
    // http://stackoverflow.com/a/1609664/148668
    if(cs == NULL)
    {
      CFStringGetCString(
        file, 
        cs_buf, 
        [(NSString*)file length]+1, 
        kCFStringEncodingUTF8);
      cs = cs_buf;
    }

    const float TRANSPARENT_BLACK[4] = {0,0,0,0};
    render_to_buffer(cs,TRANSPARENT_BLACK,bwidth,bheight,buffer);
    CGColorSpaceRef colorspace = CGColorSpaceCreateDeviceRGB();
    CFDataRef rgbData = CFDataCreate(NULL, buffer, bwidth * bheight * 4);
    CGDataProviderRef provider = CGDataProviderCreateWithCFData(rgbData);
    // http://forum.sparrow-framework.org/topic/create-uiimage-from-pixel-data-problems
    CGImageRef image = 
      CGImageCreate(
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

    if (QLPreviewRequestIsCancelled(preview))
        return noErr;

    // Draw onto context as textured rectangle
    CGContextRef cgContext = QLPreviewRequestCreateContext(preview, *(CGSize *)&size, true, NULL);
    CGRect rect = CGRectMake(0,0, width, height);
    CGContextDrawImage(cgContext, rect, image);
    QLPreviewRequestFlushContext(preview, cgContext);

    CFRelease(cgContext);
    CGImageRelease(image);
    [pool release];
    return noErr;
}


void CancelPreviewGeneration(void* thisInterface, QLPreviewRequestRef preview)
{
    // implement only if supported
}

