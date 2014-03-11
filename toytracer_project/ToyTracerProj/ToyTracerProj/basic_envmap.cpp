/***************************************************************************
* basic_envmap.cpp    (environment map plugin)                             *
*                                                                          *
* An "environment map" is a function that associates a color with any      *
* given ray, often by table-lookup (i.e. by indexing colors in a image     *
* or a high-dynamic-range map of an environment.  This "basic" version     *
* simply returns a single color, with no lookup or other computation.      *
* Its fixed color is specified as a parameter to the envmap.               *
*                                                                          *
* History:                                                                 *
*   09/27/2005  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#include "toytracer.h"
#include "params.h"

struct basic_envmap : public Envmap {
    basic_envmap() {}
    basic_envmap( const Color &c ) { color = c; }
   ~basic_envmap() {}
    virtual Color Shade( const Ray & ) const { return color; }
    virtual Plugin *ReadString( const string &params );
    virtual string MyName() const { return "basic_envmap"; }
    virtual bool Default() const { return true; }
    Color color;
    };

Plugin *basic_envmap::ReadString( const string &params ) 
    {
    Color color;
    ParamReader p( params );
    if( p["envmap"] && p[MyName()] && p[color] )
        return new basic_envmap( color );
    return NULL;
    }

REGISTER_PLUGIN( basic_envmap );



