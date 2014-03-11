/***************************************************************************
* main.h                                                                   *
*                                                                          *
* This is the main program for the "toytracer".  It reads in the camera    *
* and scene from a file, then invokes the rasterizer, which casts all the  *
* initial rays and writes the resulting image to a file.                   *
*                                                                          *
* History:                                                                 *
*   10/04/2005  Updated for 2005 graphics class.                           *
*   10/10/2004  Print registered objects, get optional file name from argv.*
*   04/03/2003  Main program now defines scene geometry.                   *
*   04/01/2003  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#include "toytracer.h"

static const string DefaultScene = "scenes/scene1";

int main( int argc, char *argv[] )
    {
    Scene  scene;
	Scene scene2;
    Camera camera;

	//whether or not to do motion blur
	bool doMotionBlur = false;

	//if motion blur is not done, this toggles whether or not to render scene 2 instead
	bool renderScene2 = false; 

	string fname = "scenes/scene1";
	string fname2 = "scenes/scene2";
	string output_fname = "scenes/scene_1_2_blur";

    // Print out a banner with the current version number of the software.

    cout << "**** ToyTracer Version "
         << TOYTRACER_MAJOR_VERSION << "." << TOYTRACER_MINOR_VERSION
         << " ****" << endl;

    // Print out the names of all the types of objects that are registered before
    // reading in the sdf file.  This list should be the first thing to check if
    // the file is not parsed correctly.

    PrintRegisteredPlugins( cout );

    // There must be at least one builder plugin to build the scene (usually by reading
    // it from a file).  Find the first that is not a "default" builder, if there is one.
    // Otherwise, use the default builder.

    Builder *builder = (Builder *)LookupPlugin( builder_plugin );
    if( builder == NULL )
        {
        cerr << "No builder plugin was found." << endl;
        return error_no_builder;
        }

    // Set the file name to a default, but substitute the argument given on the
	// command line if there is one.
    //tring fname = DefaultScene;
    //if( argc > 1 ) fname = argv[1];
	
    // Invoke the builder to construct the scene.

	if(doMotionBlur || !renderScene2){
		if( !builder->BuildScene( fname, camera, scene ) )
        {
			cerr << "Error encountered while building scene." << endl;
			return error_building_scene;
        }
	}
	if(doMotionBlur || renderScene2){
		if( !builder->BuildScene( fname2, camera, scene2 ) )
        {
			cerr << "Error encountered while building scene." << endl;
			return error_building_scene;
        }
	}

    // If a rasterizer was not specified by the builder, look to see if one has
    // been registered.

    if( scene.rasterize == NULL )
        {
        scene.rasterize = (Rasterizer *)LookupPlugin( rasterizer_plugin );
        if( scene.rasterize == NULL )
            {  
            cerr << "No rasterizer plugin was found." << endl;
            return error_no_rasterizer;
            }
        }

    // Generate the image using the rasterizer that was either specified by
    // the builder or supplied by default.

	if(doMotionBlur){
		if( !scene.rasterize->Rasterize( output_fname, camera, scene, scene2,true ) )
		{
			cerr << "Error encountered while rasterizing." << endl;
			return error_rasterizing_image;
		}
	}else if(renderScene2){
		if( !scene.rasterize->Rasterize( fname2, camera, scene2, scene,false ) )
		{
			cerr << "Error encountered while rasterizing." << endl;
			return error_rasterizing_image;
		}
	}else{
		if( !scene.rasterize->Rasterize( fname, camera, scene, scene2,false ) )
		{
			cerr << "Error encountered while rasterizing." << endl;
			return error_rasterizing_image;
		}
	}

    DestroyRegisteredPlugins();
    return no_errors;
    }

