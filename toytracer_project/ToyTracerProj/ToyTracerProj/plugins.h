/***************************************************************************
* plugins.h                                                                *
*                                                                          *
* This file defines the classes and functions that allow plugins to be     *
* created, registered, accessed, and destroyed.                            *
*                                                                          *
* History:                                                                 *
*   04/23/2003  Split off from reader.                                     *
*                                                                          *
***************************************************************************/
#ifndef __PLUGINS_INCLUDED__
#define __PLUGINS_INCLUDED__

// Here are all the plug-in types that are currently supported.
enum plugin_type {       
    primitive_plugin  = 1,  // Primitive objects, like sphere, block, cone, etc.
    aggregate_plugin  = 2,  // Aggregate objects, like list, BVH, Octree, etc.
    shader_plugin     = 3,  // A surface shader; one associated with each object.
    envmap_plugin     = 4,  // Environment maps defining what surrounds the scene.
    rasterizer_plugin = 5,  // Raserizers, which create images by tracing rays.
    builder_plugin    = 6,  // Scene builder, reads and/or constructs the scene.
    data_plugin       = 7   // A container for arbitrary data.
    };

// This is the base class of all plugins.  Each object, shader, rasterizer, etc. must
// be a subclass of this class, and they must fill in the pure virtual functions:
// ReadString, MyName, and PluginType.
struct Plugin {  
    Plugin() {}
    virtual ~Plugin() {}
    virtual Plugin *ReadString( const string & ) = NULL; // Instance from string.
    virtual string MyName() const = 0;
    virtual plugin_type PluginType() const = 0;
    virtual bool Default() const { return false; } // Use default plugins if no other defined.
    virtual void AddData( Plugin * ) {} // Optional way to pass data to a plugin.
    };

// RegisterPlugin is called by the REGISTER_PLUGIN macro.  This should never
// be called except via the macro.
extern bool RegisterPlugin(
    Plugin *plugin
    );

// PrintRegisteredPlugins writes a list of all plugins that have been reistered +
// their type to the specified output stream.  This is useful at startup, to verify 
// that all the desired plugins are present.
extern void PrintRegisteredPlugins(
    ostream &out
    );

// Return a pointer to the first plugin in the list of registered plugins
// that is of the specified type.  If "after" is set to a plugin pointer,
// this function will return a pointer to the first plugin of the specified
// type that occurs *after* that one in the list of registered plugins.  This
// allows all plugins of a given type to be accessed.
extern Plugin *LookupPlugin(
    plugin_type type,
    const Plugin *after = NULL
    );

// Call the readers of all the plugins to see if any of them recognize
// this string.  If so, create an instance of the plugin that recognizes
// the string as its parameters, and return a pointer to the newly-created
// instance as the function value.
extern Plugin *Instance_of_Plugin(
    const char *buff
    );

// DestroyRegisteredPlugins calls the destructor for each plugin.
extern void DestroyRegisteredPlugins(
    );

// The following macro is used for registering new plugins in the toytracer.
// By virtue of this macro, plugins can be added to the toytracer without modifying
// any of the existing code.  Simply define a new plugin (e.g. object, aggregate,
// shader, etc.) in a separate source file and invoke this macro using the class name
// of the new plugin.  Linking the new module into the toytracer will automatically
// add the new plugin into the list of "registered" plugins at run time.  (This macro
// uses the fact that global variables are initialized before "main" is invoked.)

#define REGISTER_PLUGIN( class_name ) \
    static bool dummy_variable_##class_name = RegisterPlugin( new class_name );

#endif


