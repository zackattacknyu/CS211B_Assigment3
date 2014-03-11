/***************************************************************************
* base.h                                                                   *
*                                                                          *
* This header file is included by every module of the system.  It is       *
* responsible for including most of the system header files needed, and    *
* for defining some fundamental structures and constants.                  *
*                                                                          *
* History:                                                                 *
*   12/11/2004  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#ifndef __BASE_INCLUDED__    // Include this file only once.
#define __BASE_INCLUDED__

// This is where (most of) the system header files are included.
// These may need to be modified on some systems.  Not all modules will
// actually need all of these headers, but it's convenient to have them
// all handled in one place.

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <cassert>

// The following names from the "std" namespace are so commonly used
// that it's convenient to refer to them without the "std:" prefix.

using std::cout;
using std::cerr;
using std::endl; 
using std::ostream;
using std::string;
using std::vector;

// We declare all the main structures here so that we can include
// pointers to these objects before their full definitions are provided.

struct Plugin;     // Things that can be added to the toytracer.
struct Object;     // Anything that can be ray traced.
struct Primitive;  // Simple geometrical objects: triangle, sphere, etc.
struct Aggregate;  // Collections of primitive and/or aggregate objects.
struct Shader;     // Shaders that can be associated with surfaces.
struct Envmap;     // An environment map, or the light "from infinity".
struct Material;   // Surface parameters available to shaders.
struct Camera;     // The parameters of the camera, image resolution, etc.
struct Scene;      // The camera, lights, object(s), etc.
struct Rasterizer; // The function that casts primary rays & creates an image.
struct Builder;    // Builds the scene, usually by reading a file (e.g. sdf).

// Miscellaneous numerical constants.

static const double 
    Pi          = 3.14159265358979,
    TwoPi       = 2.0 * Pi,
    FourPi      = 4.0 * Pi,
    DegToRad    = Pi / 180.0,  // Convert degrees to radians.
    Infinity    = 1.0E20,      // Should suffice for "infinity" in most cases.
    MachEps     = 1.0E-5,      // A reasonable value for the "machine epsilon".
    OnePlusEps  = 1.0 + MachEps,
    OneMinusEps = 1.0 - MachEps;

// Miscellaneous default values.

static const int
    default_image_width  = 400,  // Default image width (x resolution).
    default_image_height = 400,  // Default image height (y resolution).
    default_max_tree_depth = 4;  // Default cap on ray tree depth.

enum toytracer_error {
    no_errors = 0,
    error_opening_input_file,
    error_reading_input_file,
    error_opening_image_file,
    error_no_builder,
    error_building_scene,
    error_no_rasterizer,
    error_rasterizing_image
    };

#endif

