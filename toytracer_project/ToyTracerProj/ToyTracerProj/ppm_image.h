/***************************************************************************
* ppm_image.h                                                              *
*                                                                          *
* The image class defines a trivial encoding for images known as PPM       *
* format; it simply consists of an array or RGB triples, with one byte     *
* per component, preceeded by a simple header giving the size of the       *
* image.                                                                   *
*                                                                          *
* History:                                                                 *
*   10/21/2005  File i/o now based on ofstream.                            *
*   04/01/2003  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#ifndef __PPM_IMAGE_INCLUDED__
#define __PPM_IMAGE_INCLUDED__

#include <string>

typedef unsigned char channel;  

struct Pixel {
    Pixel() { r = 0; g = 0; b = 0; }
    Pixel( channel R, channel G, channel B ) { r = R; g = G; b = B; }
    channel r;
    channel g;
    channel b;
    };

struct PPM_Image {
    PPM_Image( int x_res, int y_res );
   ~PPM_Image() { delete[] pixels; }
    bool Read ( std::string file_name );
    bool Write( std::string file_name );
    inline  Pixel &operator()( int i, int j ) { return *( pixels + ( i * width + j ) ); }  
    Pixel  *pixels;
    int     width;
    int     height;
    };

extern bool Overwrite_PPM_Image(
    std::string file_name
    );

#endif


