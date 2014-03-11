/***************************************************************************
* basic_rasterizer.cpp    (rasterizer plugin)                              *
*                                                                          *
* A "rasterizer" is responsible for tracing all the primary rays needed to *
* create an image, filling in the resulting matrix of color values (i.e.   *
* the "raster"), and saving the results as an image file.  This "basic"    *
* rasterizer simply casts a single ray per pixel and saves the result as   *
* a PPM image.  It does no anti-aliasing.                                  *
*                                                                          *
* History:                                                                 *
*   10/03/2005  Made rasterizer a plugin.  Line numbers written in place.  *
*   12/19/2004  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#include "toytracer.h"
#include "ppm_image.h"
#include "params.h"
#include "util.h"

/*
*Deciding on how many rays is tricky. Too few and there is not much anti-aliasing.
*	Too many and it takes too long to compute. 
*2000 is a good number once you can wait a while. It should be lowered to 20 though
*	when you are testing the code. 
*/
static const double numRaysAntiAliasing = 1;
static const double numRaysDepthOfField = 1;

struct basic_rasterizer : public Rasterizer {
    basic_rasterizer() {}
    virtual ~basic_rasterizer() {}
    virtual bool Rasterize( string fname, const Camera &, const Scene & , const Scene &, const bool &doMotionBlur) const;
    virtual Plugin *ReadString( const string &params );
    virtual string MyName() const { return "basic_rasterizer"; }
    virtual bool Default() const { return true; }
    };

REGISTER_PLUGIN( basic_rasterizer );

Plugin *basic_rasterizer::ReadString( const string &params ) 
    {
    ParamReader p( params );
    if( p["rasterizer"] && p[MyName()] ) return new basic_rasterizer();
    return NULL;
    }

// This is a trivial tone mapper; it merely maps values that are
// in [0,1] and maps them to integers between 0 and 255.  If the
// real value is above 1, it merely truncates.  A true tone mapper
// would attempt to handle very large values nicely, without
// truncation; that is, it would try to compensate for the fact that
// displays have a very limited dynamic range.
static Pixel ToneMap( const Color &color )
    {
    int red   = (int)floor( 256 * color.red   );
    int green = (int)floor( 256 * color.green );
    int blue  = (int)floor( 256 * color.blue  );
    channel r = (channel)( red   >= 255 ? 255 : red   ); 
    channel g = (channel)( green >= 255 ? 255 : green ); 
    channel b = (channel)( blue  >= 255 ? 255 : blue  );
    return Pixel( r, g, b );
    }

// Rasterize casts all the initial rays starting from the eye.
// This trivial version simply casts one ray per pixel, in raster
// order, then writes the pixels out to a file.
bool basic_rasterizer::Rasterize( string file_name, const Camera &cam, const Scene &scene, const Scene &scene2, const bool &doMotionBlur ) const
    {
    file_name += ".ppm";

    // Make sure the file is accessible by overwriting it now.  That way, if the file
    // is not accessible, we'll find out now instead of waiting until the image is ready
    // to be written.

    if( !Overwrite_PPM_Image( file_name ) )
        {
        cerr << "Error: Could not open file " << file_name << " for writing." << endl;
        return false;
        }

    // Create an image of the given resolution.

    PPM_Image I( cam.x_res, cam.y_res );

    // Initialize all the fields of the first-generation ray except for "direction".

    Ray ray;
    ray.origin     = cam.eye;     // All initial rays originate from the eye.
    ray.type       = generic_ray; // These rays are given no special meaning.
    ray.generation = 1;           // Rays cast from the eye are first-generation.

    const double xmin   = cam.x_win.min;
    const double ymax   = cam.y_win.max;
    const double width  = Len( cam.x_win );
    const double height = Len( cam.y_win );

    // Compute increments etc. based on the camera geometry.  These will be used
    // to define the ray direction at each pixel.

    const Vec3 G ( Unit( cam.lookat - cam.eye ) );          // Gaze direction.
    const Vec3 U ( Unit( cam.up / G ) );                    // Up vector.
    const Vec3 R ( Unit( G ^ U ) );                         // Right vector.
    const Vec3 O ( cam.vpdist * G + xmin * R + ymax * U );  // "Origin" of the 3D raster.
    const Vec3 dR( width  * R / cam.x_res );                // Right increments.
    const Vec3 dU( height * U / cam.y_res );                // Up increments.

    // Loop over the entire image, casting a single ray per pixel.  This nested loop
    // must be modified to accommodate anti-aliasing.
	
	Color currentColor = Color();
	double randomX = 0.5;
	double randomY = 0.5;
	double randomDU = 0.5;
	double randomDR = 0.5;
	const double focalLength = 5;
	Vec3 imagePlanePoint;
	Vec3 currentOrigin;
	Vec3 currentDirection;
	const double radius = 20;
	
    cout << "Rendering line 0";
    for( unsigned int i = 0; i < cam.y_res; i++ )
        {
        // Overwrite the line number written to the console.
        cout << rubout( i ) << (i+1);
        cout.flush();
        for( unsigned int j = 0; j < cam.x_res; j++ )
            {
				currentColor = Color();

				//shoots multiple rays in the pixel window
				for(int rayNum = 0; rayNum < numRaysAntiAliasing; rayNum++){

					//generates a random pair in [0,1]x[0,1] to be used as the current ray
					if(numRaysAntiAliasing > 1){ //in case we are not doing anti-aliasing
						randomX = (double)rand() / RAND_MAX;
						randomY = (double)rand() / RAND_MAX;
					}
					

					//shoots the random ray found and gets its color
					ray.direction = Unit( O + (j + randomX) * dR - (i + randomY) * dU  );
					imagePlanePoint = cam.eye + focalLength*ray.direction;

					//shoot the ray from different origin points for depth of field effect
					for(int dofNum = 0; dofNum < numRaysDepthOfField; dofNum++){

						if(numRaysDepthOfField > 1){ //in case there is no depth of field
							randomDU = (double)rand() / RAND_MAX;
							randomDR = (double)rand() / RAND_MAX;
						}
						
						//jitter the camera position
						ray.origin = cam.eye + (randomDU-0.5)*radius*dR + (randomDR-0.5)*radius*dU;
					
						//get the direction from the jittered position to the image plane position
						ray.direction = Unit(imagePlanePoint - ray.origin);

						
						if(doMotionBlur){
							currentColor = currentColor + 0.15*scene.Trace(ray) + 0.85*scene2.Trace(ray);
						}else{
							currentColor = currentColor + scene.Trace(ray);
						}

					}
					
				}

				//blends the colors together of the found rays
				currentColor = currentColor/(numRaysAntiAliasing*numRaysDepthOfField);

				I(i,j) = ToneMap(currentColor);
            }
        }

    // Thus far the image exists only in memory.  Now write it out to a file.

    cout << "\nWriting image file " << file_name << "... ";
    cout.flush();
    I.Write( file_name );
    cout << "done." << endl;
    return true;
    }