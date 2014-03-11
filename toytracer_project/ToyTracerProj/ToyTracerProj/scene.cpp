/***************************************************************************
* scene.cpp                                                                *
*                                                                          *
* The "scene" structure defines two of the most fundamental functions of   *
* the ray tracer: "Cast" and "Trace".  Both are rather trivial, as the     *
* real computations are encapsulated in the actual objects comprising      *
* the sceen.                                                               *
*                                                                          *
* History:                                                                 *
*   09/29/2005  Updated for 2005 graphics class.                           *
*   10/16/2004  Check for "ignored" object in "Cast".                      *
*   09/29/2004  Updated for Fall 2004 class.                               *
*   04/14/2003  Point lights are now point objects with emission.          *
*   04/01/2003  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#include "toytracer.h"

static const Color
    default_background_color  = Color( 0.15, 0.25, 0.35 ),
    default_termination_color = Color( 0.15, 0.15, 0.15 );

Scene::Scene()
    {
    object    = NULL;
    envmap    = NULL; 
    rasterize = NULL;
    max_tree_depth = default_max_tree_depth;
    }

// Cast finds the first point of intersection (if there is one)
// between a ray and a list of geometric objects.  If no intersection
// exists, the function returns false.  Information about the
// closest object hit is returned in "hitinfo".  Note that the single
// object associated with the scene will typically be an aggregate (so
// that the scene can contain many primitive objects).

bool Scene::Cast( const Ray &ray, HitInfo &hitinfo ) const
    {
    if( object == NULL || object == hitinfo.ignore ) return false;
    if( object->Intersect( ray, hitinfo ) )
        {
        hitinfo.ray = ray; // Save the ray in world coordinates.
        return true;
	    }
    return false;
    }

// Trace is the most fundamental of all the ray tracing functions.  It
// answers the query "What color do I see looking along the given ray
// in the current scene?"  This is an inherently recursive process, as
// trace may again be called as a result of the ray hitting a reflecting
// object.  To prevent the possibility of infinite recursion, a maximum
// depth is placed on the resulting ray tree.

Color Scene::Trace( const Ray &ray ) const
    {
    Color   color;               // The color to return.
    HitInfo hitinfo;             // Holds info to pass to shader.
    hitinfo.ignore = NULL;       // Don't ignore any objects.
    hitinfo.distance = Infinity; // Follow the full ray.

    if( ray.generation > max_tree_depth )
        {
        // The recursion has bottomed out.  Do something!
        color = default_termination_color;
        }
    else if( Cast( ray, hitinfo ) )
        {
        if( hitinfo.object == NULL ) return Green;
        // The ray hits an object, so shade the point that the ray hit.
        // Cast has put all necessary information for Shade into "hitinfo".
        // Use the shader associated with the object, if there is one.
        Shader *shader = hitinfo.object->shader;
        if( shader != NULL )
             color = shader->Shade( *this, hitinfo );   // Use the associated shader.
        else color = hitinfo.object->material->diffuse; // Use the diffuse color.
        }
    else 
        {
        // The ray has failed to hit anything.  Use the environment map
        // associated with the scene to determine the color (if there is one).
        // If no envmap was specified, use the default background color.
        const Envmap *env = envmap;
        if( ray.from != NULL && ray.from->envmap != NULL ) env = ray.from->envmap;
        if( env != NULL ) 
             color = env->Shade( hitinfo.ray ); // Use the associated env.
        else color = default_background_color;  // Use the default color.
        }
    
    return color;
    }

