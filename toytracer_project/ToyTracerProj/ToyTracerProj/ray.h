/***************************************************************************
* ray.h                                                                    *
*                                                                          *
* The "Ray" structure defines a ray in 3-space, consisting of an origin    *
* (denoted by "Q") and a direction (denoted by "R").  Several additional   *
* fields are also included (e.g. type and generation) as a convenience for *
* some ray tracing algorithms.                                             *
*                                                                          *                                                                        
* History:                                                                 *
*   12/11/2004  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#ifndef __RAY_INCLUDED__
#define __RAY_INCLUDED__

#include "base.h"
#include "vec3.h"

// Possible ray "types" that may be set and used by a shader and/or an
// acceleration method.
enum ray_type {         // A flag that may affect the processing of a ray.
    undefined_ray = 0,  // Nothing known about this ray.
    primary_ray   = 1,  // A ray cast from the eye by the rasterizer.
    generic_ray   = 2,  // Any ray; No special meaning.
    shadow_ray    = 3,  // A ray cast for shadow detection only.
    indirect_ray  = 4,  // A ray cast from a surface to sample illumination.
    light_ray     = 5,  // A ray cast from a light source (e.g. for photon mapping).
    special_ray   = 6   // A ray used for some other special purpose.
    };

struct Ray {            // A ray in R3.
    inline Ray();
    inline Ray( const Ray &r );
    virtual ~Ray() {}
    Vec3 origin;         // The ray originates from this point.
    Vec3 direction;      // Unit vector indicating direction of ray.
    int  type;           // Different rays may be processed differently.
	double ref_index;	//ref index of current material the ray is inside
    unsigned generation; // How deep in the ray tree.  1 == generated from eye.
    const Object *from;  // The object from which the ray was cast.
    };

inline Ray::Ray()
    {
    generation = 1;
    type = generic_ray;
    from = NULL;
	ref_index = 1.0;
    }

inline Ray::Ray( const Ray &r )
    {
    origin     = r.origin;
    direction  = r.direction;
    generation = r.generation;
    type       = r.type;;
    from       = r.from;
	ref_index = 1.0;
    }


// Compute the reflected ray given the incident ray (i.e. directed
// toward the surface), and the normal to the surface.  The normal
// may be directed away from or into the surface.  Both the surface
// normal and the ray direction vector are assumed to be normalized.
inline Vec3 Reflect( const Ray &r, const Vec3 &N )
    {
    Vec3 U( r.direction );
    return U - ( 2.0 * ( U * N ) ) * N;
    }

// Compute the refracted ray given the incident ray (i.e. directed
// toward the surface), the normal to the surface, and the ratio of the
// refractive indices of the material containing the ray and the
// material into which the ray is refracted.  The normal
// may be directed away from or into the surface.  Both the surface
// normal and the ray direction vector are assumed to be normalized.
inline Vec3 Refract( const Ray &r, const Vec3 &N, double eta1_over_eta2 )
    {
    Vec3 U( r.direction );
    return U - ( 2.0 * ( U * N ) ) * N;
    }

// An output method, useful for debugging.
inline ostream &operator<<( ostream &out, const Ray &r )
    {
    out << "[ray:"
        << " org "  << r.origin
        << " dir "  << r.direction
        << " type " << r.type
        << " gen "  << r.generation
        << "]";
    return out;
    };

#endif

