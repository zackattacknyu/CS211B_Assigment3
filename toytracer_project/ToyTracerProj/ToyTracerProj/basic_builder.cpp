/***************************************************************************
* basic_builder.cpp                                                        *
*                                                                          *
* This file defines a basic BuildScene function, which reads a simple text *
* description of a scene and the camera.                                   *
*                                                                          *
* History:                                                                 *
*   04/23/2006  The reader is now a "Builder" plugin.                      *
*   09/29/2005  Updated for 2005 class.                                    *
*   10/23/2004  Changed handling of default colors.                        *
*   10/16/2004  Minor updates (e.g. added more error messages).            *
*   10/06/2004  Added 3x4 matrix reader & support for aggregate objects.   *
*   04/04/2003  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#include <fstream>
#include "toytracer.h"
#include "util.h"
#include "params.h"

struct basic_builder : public Builder {
    basic_builder() {}
    virtual ~basic_builder() {}
    virtual bool BuildScene( string command, Camera &, Scene & ) const;
    virtual Plugin *ReadString( const string &params );
    virtual string MyName() const { return "basic_builder"; }
    virtual bool Default() const { return true; }
    };

REGISTER_PLUGIN( basic_builder );

Plugin *basic_builder::ReadString( const string &params ) 
    {
    ParamReader p( params );
    if( p["builder"] && p[MyName()] ) return new basic_builder();
    return NULL;
    }

static bool Skip( char *line )
    {
    for( char c = *line; c != '\0' && c != '\n' && c != '#'; c = *++line )
        {
        if( c != ' ' && c != '\t' ) return false;
        }
    return true;
    }

static Material *Copy( Material* &mat, const Material &material )
    {
    if( mat == NULL || !(*mat == material) ) mat = new Material( material );
    return mat;
    }

// This is a very minimal scene description reader.  It assumes that
// each line contains a complete entity: an object definition, or
// a camera parameter, or a material parameter, etc.  (Blank lines, and
// lines that begin with "#" are also allowed.)  It creates object instances,
// links them together, and fills in the fields of the scene and camera as it
// reads successive lines of the file.  Lines may be padded with blanks or
// tabs to increase readability.

bool basic_builder::BuildScene( string file_name, Camera &camera, Scene &scene ) const
    {
    int line_num = 0;
    char input_line[512];
    Plugin    *dat = NULL;  // A container for arbitrary data.
    Object    *obj = NULL;  // The current object, which was the last object created.
    Shader    *shd = NULL;  // The current shader. 
    Aggregate *agg = NULL;  // Current aggregate object, to which all objects are now added.
    Envmap    *env = NULL;  // Current environment map.
    Material  *mat = NULL;  // Current material pointer.
    Material   material;    // Current material.

    scene.object    = NULL;
    scene.envmap    = NULL;
    scene.rasterize = NULL;

    // Attempt to open the input file.

    file_name += ".sdf";
    std::ifstream fin;
    fin.open( file_name.c_str() );
    if( fin.fail() )
        {
        cerr << "Error: Could not open file " << file_name << endl;
        return false;  // Report failure.
        }
    cout << "Reading " << file_name << "... ";
    cout.flush();

    // Set some defaults.

    material.diffuse      = White;
    material.emission     = Black;
    material.specular     = White;
    material.ambient      = Black;
    material.reflectivity = Black;
    material.translucency = Black;
    material.ref_index    = 0.0;
    material.Phong_exp    = 0.0;

    camera.x_res = default_image_width;
    camera.y_res = default_image_height;
    camera.x_win = Interval( -1.0, 1.0 );
    camera.y_win = Interval( -1.0, 1.0 );

    // Process lines until the end of file is reached.
    // Print a warning for all lines that are unrecognizable.

    while( fin.getline( input_line, 512 ) ) 
        {
        line_num++;
        if( Skip( input_line ) ) continue;

        // Ask each registered object if it recognizes the line.  If it does, it will
        // create a new instance of the object and return it as the function value.

        Plugin *plg = Instance_of_Plugin( input_line );

        if( plg != NULL )
            {
            switch( plg->PluginType() )
                {
                case data_plugin:
                    if( obj == NULL ) cerr << "Error: data ignored.  Line " << line_num << endl;
                    else obj->AddData( obj );
                    break;

                case shader_plugin:
                    shd = (Shader*)plg;
                    break;					
 
                case aggregate_plugin:
                    obj = (Object*)plg;
                    obj->shader   = shd;
                    obj->envmap   = env;
                    obj->material = Copy( mat, material );
                    obj->parent   = agg;
                    if( Emitter( material ) )
                        {
                        cerr << "Error: An aggregate object cannot be an emitter.  Line "
                             << line_num << ": "
                             << input_line << endl;
                        return false;
                        }
                    if( agg != NULL ) // If there is alrealy an agg obj, this one is a child.
                        {
                        // agg->AddChild( obj );
                        // Do not add aggregates as children until they are complete.
                        agg->material = Copy( mat, material );
                        }
                    else if( scene.object == NULL ) scene.object = obj;
	            agg = (Aggregate *)obj;
                    break;

                case primitive_plugin:
                    obj = (Object*)plg;
                    obj->shader   = shd;
                    obj->envmap   = env;
                    obj->material = Copy( mat, material );
                    obj->parent   = agg;
                    if( Emitter( material ) ) scene.lights.push_back( obj );
                    if( agg != NULL )
                        {
                        agg->AddChild( obj );
                        agg->material = Copy( mat, material );
                        }
                    else if( scene.object == NULL ) scene.object = obj;
                    break;

                case envmap_plugin:
                    env = (Envmap*)plg;
                    // If an environment map is set before any object is created, use it as
                    // the background.
                    if( obj == NULL ) scene.envmap = env;
                    break;

                case rasterizer_plugin:
                    if( scene.rasterize != NULL )
                        {
                         cerr << "Error: More than one rasterizer specified.  Line "
                             << line_num << ": "
                             << input_line << endl;
                        return false;
                        }
                    scene.rasterize = (Rasterizer *)plg;
                    break;
                } 
            continue;
            }

        // Now look for all the other stuff...  materials, camera, lights, etc.

        ParamReader get( input_line );
 
        if( get["eye"]          && get[camera.eye]            ) continue;
        if( get["lookat"]       && get[camera.lookat]         ) continue;            
        if( get["up"]           && get[camera.up]             ) continue;            
        if( get["vpdist"]       && get[camera.vpdist]         ) continue;                   
        if( get["x_res"]        && get[camera.x_res]          ) continue;
        if( get["y_res"]        && get[camera.y_res]          ) continue;
        if( get["x_win"]        && get[camera.x_win]          ) continue;
        if( get["y_win"]        && get[camera.y_win]          ) continue;
        if( get["ambient"]      && get[material.ambient]      ) continue;            
        if( get["diffuse"]      && get[material.diffuse]      ) continue;
        if( get["specular"]     && get[material.specular]     ) continue;
        if( get["emission"]     && get[material.emission]     ) continue;
        if( get["reflectivity"] && get[material.reflectivity] ) continue;
        if( get["translucency"] && get[material.translucency] ) continue;
        if( get["Phong_exp"]    && get[material.Phong_exp]    ) continue;
        if( get["ref_index"]    && get[material.ref_index]    ) continue;

        // If no object is defined at this point, it's an error.

        if( obj == NULL )
            {
            cerr << "Error reading scene file; No object defined at line "
                 << line_num << ": "
                 << input_line << endl;
            return false;
            }

        // Look for an end statement that closes the current aggregate.  This allows us to nest aggregates.

        if( get["end"] )
            { 
            if( agg == NULL )
                {
                cerr << "Error: end statement found outside an aggregate object at line "
                     << line_num << ": "
                     << input_line << endl;
                return false;
                }
            
            // Go back to adding objects to the parent object (if there is one).

            agg->Close(); // Signal the aggregate that it is now complete.
            Object *closed_agg = agg;
            agg = agg->parent;
            if( agg != NULL )
                {
                material = *(agg->material);
                mat = agg->material;
                shd = agg->shader;
                env = agg->envmap;
                agg->AddChild( closed_agg ); // Add the agg that just ended.
                }
            continue;
            }
        
        // If nothing matched, it's an error.  Print a warning and continue processing.

        cerr << "Warning: Unrecognized or ill-formed command at line "
             << line_num << ": "
             << input_line << endl;
        }

    if( agg != NULL )
        {
        cerr << "Error: Top-most aggregate object has no 'end' statement." << endl;
        return false;
        }

    cout << "done." << endl;
    return true;
    }
