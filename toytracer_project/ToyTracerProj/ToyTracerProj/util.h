/***************************************************************************
* util.h                                                                   *
*                                                                          *
* Miscellaneous utilities, such as predicates on materials & objects.      *
*                                                                          *
* History:                                                                 *
*   12/11/2004  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#ifndef __UTIL_INCLUDED__
#define __UTIL_INCLUDED__

#include "toytracer.h"

inline double sqr( double x )
    {
    return x * x;
    }

inline double min( double x, double y )
    {
    return x <= y ? x : y;
    }

inline double max( double x, double y )
    {
    return x >= y ? x : y;
    }

inline bool Emitter( const Object *obj )
    { 
    return ( obj != NULL ) & ( obj->material->emission != 0.0 );
    }

inline bool Emitter( const Material &mat )
    { 
    return mat.emission != 0.0;
    }

inline bool Reflective( const Object *obj )
    { 
    return ( obj != NULL ) & ( obj->material->reflectivity != 0.0 );
    }

inline bool Translucent( const Object *obj )
    { 
    return ( obj != NULL ) & ( obj->material->translucency != 0.0 );
    }

extern const char *rubout(  // Return enough backspaces to rub out the integer i.
    int i
    );

extern AABB GetBox(  // Construct a box from three slabs.
    const Object &obj
    );

extern bool operator==(
    const Material &a,
    const Material &b
    );

extern double min(
    double x,
    double y,
    double z
    );

extern double max(
    double x,
    double y,
    double z
    );

extern double rand(    // Return a random number uniformly distributed in [a,b].
    double a,
    double b
    );

extern string ToString(
    plugin_type ptype
    );

#endif

