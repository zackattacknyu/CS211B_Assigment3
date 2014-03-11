/***************************************************************************
* toytracer.h                                                              *
*                                                                          *
* This is the header file for a "Toy" ray tracer, or "toytracer".  Every   *
* module of the system includes this header as it defines all of the       *
* fundamental structures needed by the ray tracer.                         *
*                                                                          *
* History:                                                                 *
*   10/18/2005  Added Item & Primitive base classes.
*   09/29/2005  Now supports more plugins, including shaders.              *
*   10/10/2004  Added Aggregate sub-class & REGISTER_OBJECT macro.         *
*   10/06/2004  Added type to ray & ray to HitInfo.  Removed HitGeom.      *
*   09/29/2004  Updated for Fall 2004 class.                               *
*   04/10/2003  Added GetSamples to Object class.                          *
*   04/06/2003  Added Object class.                                        *
*   04/01/2003  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#ifndef __TOYTRACER_INCLUDED__  // Include this file only once.
#define __TOYTRACER_INCLUDED__

#define TOYTRACER_MAJOR_VERSION 6  // Indicates incompatible versions.
#define TOYTRACER_MINOR_VERSION 0  // Indicates minor fixes & extensions.

#include "base.h"         // Includes system headers and defines constants.
#include "vec3.h"         // Defines the Vec3 class, which are points in R3.
#include "vec2.h"         // Defines the Vec2 class, which are points in R2.
#include "mat3x3.h"       // Defines 3x3 matrices.
#include "mat3x4.h"       // Defines 3x4 matrices; i.e. affine transforms.
#include "color.h"        // Defines the Color class; real RGB values.
#include "interval.h"     // Defines a (min,max) interval of the real line.
#include "aabb.h"         // Defines a 3D axis-aligned bounding box.
#include "ray.h"          // Defines rays in 3-space: origin, direction, etc.
#include "plugins.h"      // Defines functions for accessing plugins.

struct Material {         // Surface material for shading.
    Color  diffuse;       // Diffuse color.
    Color  specular;      // Color of highlights.
    Color  emission;      // Emitted light.
    Color  ambient;       // Ambient light (from all directions).
    Color  reflectivity;  // Weights for refleced light, between 0 and 1.
    Color  translucency;  // Weights for refracted light, between 0 and 1.
    double Phong_exp;     // Phong exponent for specular highlights.
    double ref_index;     // Refractive index.
    long   type;          // Reserved for future use.
    };

struct HitInfo {          // Records all info at ray-object intersection.
    const Object *ignore; // One object in scene to ignore (used by Cast).
    const Object *object; // The object that was hit (set by Intersect).
    double  distance;     // Distance to hit (used & reset by Intersect).
    Vec3    point;        // ray-object intersection point (set by Intersect).
    Vec3    normal;       // Surface normal (set by Intersect).
    Vec2    uv;           // Texture coordinates (set by intersect).
    Ray     ray;          // The ray that hit the surface (set by Cast).
    };

struct Sample {           // A point and weight returned from a sampling algorithm.
    Vec3   P;
    double w;
    };

struct Camera {           // Defines the position of the eye/camera.
    Vec3     eye;         // Position of eye.
    Vec3     lookat;      // The point we are looking toward.
    Vec3     up;          // A vector not parallel to the direction of gaze.
    double   vpdist;      // Distance to the view plane.
    Interval x_win;       // Horizontal extent of view window (typically [-1,1]).
    Interval y_win;       // Vertical extent of view window (typically [-1,1]).
    unsigned x_res;       // Horizontal image resolution in pixels.
    unsigned y_res;       // Vertical image resolution in pixels.
    };

struct Scene {
    Scene();
   ~Scene() { lights.clear(); }
    Color Trace( const Ray &ray ) const;
    bool  Cast ( const Ray &ray, HitInfo &hitinfo ) const;
    virtual const Object *GetLight( unsigned i ) const { return lights[i]; } 
    virtual unsigned NumLights() const { return lights.size(); }
    Envmap     *envmap;      // Global environment map, if ray hits nothing. 
    Object     *object;      // A single primitve or an aggregate object.
    Rasterizer *rasterize;   // This casts all primary rays & makes the image.
    vector<Object*> lights;  // All objects that are emitters.  
    unsigned max_tree_depth; // Limit on depth of the ray tree.
    };

struct Shader : Plugin {  // Each object has an associated shader.
    Shader() {}
    virtual ~Shader() {}
    virtual Color Shade( const Scene &scene, const HitInfo &hitinfo ) const = 0;
    virtual plugin_type PluginType() const { return shader_plugin; }
    };

struct Envmap : Plugin { // Each object can have its own environment map.
    Envmap() {}
    virtual ~Envmap() {}
    virtual Color Shade( const Ray &ray ) const = 0;
    virtual plugin_type PluginType() const { return envmap_plugin; }
    };

struct Object : Plugin { // Base class for all objects that can be ray traced.
    Object() { material = 0; shader = 0; envmap = 0; parent = 0; }
    virtual ~Object() {}
    virtual bool Intersect( const Ray &ray, HitInfo & ) const = 0;
    virtual unsigned GetSamples( const Vec3 &P, const Vec3 &N, Sample *samples, unsigned n ) const { return 0; }
    virtual bool Inside( const Vec3 & ) const = 0;
    virtual Interval GetSlab( const Vec3 & ) const = 0;
    virtual double Cost() const { return 1.0; }
    Material  *material;
    Shader    *shader;
    Envmap    *envmap;
    Aggregate *parent; 
    };

struct Primitive : Object { // Base class for all primitive objects.
    Primitive() {}
    virtual ~Primitive() {}
    virtual plugin_type PluginType() const { return primitive_plugin; }
    };

struct Aggregate : Object {  // Base class for all Aggregate objects.
    Aggregate() {}
    virtual ~Aggregate() { children.clear(); }
    virtual void AddChild( Object *obj ) { children.push_back( obj ); }
    virtual void Close() {} // Called when all children have been added.
    virtual unsigned NumChildren() const { return children.size(); }
    virtual const Object *GetChild( unsigned i ) const { return children[i]; } 
    virtual plugin_type PluginType() const { return aggregate_plugin; }
    vector <Object*> children;
    };

struct Builder : Plugin {  // Creates the scene to be rendered.
    Builder() {}
    virtual ~Builder() {}
    virtual bool BuildScene(
		string command,   // Typically, this is the file to read.
        Camera &camera,   // Returns the view.
        Scene  &scene     // Returns scene description: object, envmap, etc.
        ) const = 0;
    virtual plugin_type PluginType() const { return builder_plugin; }
    };

struct Rasterizer : Plugin {  // The rasterizer creates all the primary rays.
    Rasterizer() {}
    virtual ~Rasterizer() {}
    virtual bool Rasterize(
        string fname,         // File to write to (must include the extension).
        const Camera &camera, // Defines the view.
        const Scene &scene,    // Global scene description: object, envmap, etc.
		const Scene &scene2,
		const bool &doMotionBlur
        ) const = 0;
    virtual plugin_type PluginType() const { return rasterizer_plugin; }
    };

#endif

