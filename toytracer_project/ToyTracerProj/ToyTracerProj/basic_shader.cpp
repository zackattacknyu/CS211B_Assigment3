/***************************************************************************
* basic_shader.cpp   (shader plugin)                                       *
*                                                                          *
* This file defines a very simple ray tracing shader for the toy tracer.   *
* The job of the shader is to determine the color of the surface, viewed   *
* from the origin of the ray that hit the surface, taking into account the *
* surface material, light sources, and other objects in the scene.         *                          *
*                                                                          *
* History:                                                                 *
*   10/03/2005  Updated for Fall 2005 class.                               *
*   09/29/2004  Updated for Fall 2004 class.                               *
*   04/14/2003  Point lights are now point objects with emission.          *
*   04/01/2003  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#include "toytracer.h"
#include "util.h"
#include "params.h"

struct basic_shader : public Shader {
    basic_shader() {}
   ~basic_shader() {}
    virtual Color Shade( const Scene &, const HitInfo & ) const;
    virtual Plugin *ReadString( const string &params );
    virtual string MyName() const { return "basic_shader"; }
    virtual bool Default() const { return true; }
	virtual Vec3 RefractionDirection(double n_1, double n_2, Vec3 incomingVector,Vec3 normalVector) const;
    };

REGISTER_PLUGIN( basic_shader );

Plugin *basic_shader::ReadString( const string &params ) 
    {
    ParamReader get( params );
    if( get["shader"] && get[MyName()] ) return new basic_shader();
    return NULL;
    }


Color basic_shader::Shade( const Scene &scene, const HitInfo &hit ) const
    {
    Ray ray;
	Ray reflectionRay;
	Ray refractedRay;
	Ray secondaryRefractedRay;

    HitInfo otherhit;
	HitInfo refractionHit;
    static const double epsilon = 1.0E-6;
    if( Emitter( hit.object ) ) return hit.object->material->emission;

    Material *mat   = hit.object->material;
    Color  diffuse  = mat->diffuse;
    Color  specular = mat->specular;
    Color  color    = mat->ambient * diffuse;
    Vec3   O = hit.ray.origin;
    Vec3   P = hit.point;
    Vec3   N = hit.normal;
    Vec3   E = Unit( O - P );
    Vec3   R = Unit( ( 2.0 * ( E * N ) ) * N - E );
    Color  r = mat->reflectivity;
    double e = mat->Phong_exp;
    double k = mat->ref_index;
	Color t = mat->translucency;

	if(hit.object->Inside(P)){
		P = P + epsilon*N;
	}

    if( E * N < 0.0 ) N = -N;  // Flip the normal if necessary.

	//get the attentuation
	//	A = 1/(a + b*r + c*r^2)
	double attenuation_a = 0.0;
	double attenuation_b = 0.0;
	double attenuation_c = 0.02;
	double attenuation;
	double lightDistance;
	Vec3 lightVector;
	double diffuseFactor;
	double specularFactor;
	Color diffuseColor = Color();
	Color specularColor = Color();
	Color finalColor;
	Color colorWithLighting;
	Color reflectedColor;
	Color refractedColor;
	Vec3 currentR;
	bool objectWasHit = false;

    for( unsigned i = 0; i < scene.NumLights(); i++ )
        {
        const Object *light = scene.GetLight(i);
        Color emission = light->material->emission;
        AABB box = GetBox( *light );
        Vec3 LightPos( Center( box ) ); 

		//gets the light Vector
		lightVector = LightPos - P;
		lightDistance = Length(lightVector);
		lightVector = Unit(lightVector);
		
		//gets the attenuation factor
		attenuation = 1/(attenuation_a + attenuation_b*lightDistance + attenuation_c*lightDistance*lightDistance);

		//light ray to case to determine occulsion
		ray.origin = P;
		ray.direction = lightVector;
		HitInfo objectHit;
		//objectHit.ignore = NULL;
		objectHit.distance = Infinity;

		const int numRaysSoftShadows = 1;
		double shadowFactor = 0;
		double randomLightDeltaY;
		double randomLightDeltaZ;
		Vec3 deltaVector;
		Vec3 currentLightVector;

		//calculate the soft shadow
		for(int rayIndex = 0; rayIndex < numRaysSoftShadows; rayIndex++){


			if(numRaysSoftShadows > 1){
				//generates two numbers between -0.05 and 0.05
				randomLightDeltaY = ( ((double)rand() / RAND_MAX) - 0.5)/10.0;
				randomLightDeltaZ = ( ((double)rand() / RAND_MAX) - 0.5)/10.0;
				deltaVector = Vec3(0.0,randomLightDeltaY,randomLightDeltaZ);
			}else{
				deltaVector = Vec3(0.0,0.0,0.0);
			}
			
			currentLightVector = lightVector + deltaVector;

			ray.direction = currentLightVector;
			//objectHit.ignore = NULL;
			objectHit.distance = Infinity;

			if(scene.Cast(ray,objectHit) ){
				if(objectHit.object != NULL){
					shadowFactor = shadowFactor + 1;
				}
			}
		}

		//calculate the shadow factor
		shadowFactor = shadowFactor/numRaysSoftShadows;
		shadowFactor = 1-shadowFactor;
		
		//gets the diffuse component
		diffuseFactor = max(0,lightVector*N);
		
		//gets the specular component
		currentR = Unit(2.0*(N*lightVector)*N - lightVector);
		specularFactor = max(0, pow(currentR*E,e) );

		if(diffuseFactor == 0){
			specularFactor = 0;
		}

		//calculate the new color
		specularColor = specularColor + shadowFactor*(attenuation*specularFactor)*emission;
		diffuseColor = diffuseColor + shadowFactor*(attenuation*diffuseFactor)*emission;
        }

	colorWithLighting = color + diffuseColor*diffuse + specularColor*specular;

	//set variables for reflection
	reflectionRay.origin = P;
	reflectionRay.direction = R;
	reflectionRay.generation = hit.ray.generation + 1;
	reflectedColor = Color();

	//set variables for refraction

	refractedRay.origin = P-epsilon*N;

	Vec3 refractionDir = RefractionDirection(1.0,k,E,N);
	refractedRay.direction = refractionDir; //for refraction
	refractedRay.generation = hit.ray.generation + 1;
	refractedColor = Color();
 	refractionHit.distance = Infinity;

	//only do refraction if the transluency is greater than zero
	//	this is an optimization so unnecessary refractions are not calculated
	if( (t.red + t.green + t.blue) > epsilon){
		if(scene.Cast(refractedRay,refractionHit)){

			bool insideMaterial = false;
			Vec3 currentNormal;
			Vec3 previousRefractedDirection;

			do{
				currentNormal = Unit(refractionHit.normal);
				previousRefractedDirection = refractedRay.direction;
				refractedRay.direction = RefractionDirection(k,1.0,-1*refractedRay.direction,-1*currentNormal);

				if(Length(refractedRay.direction) < epsilon){
					insideMaterial = true;
					refractionHit.distance = Infinity;
					refractedRay.origin = refractionHit.point - epsilon*currentNormal;
					refractedRay.direction = Unit( ( 2.0 * ( previousRefractedDirection * currentNormal ) ) * (-1*currentNormal) + previousRefractedDirection );


					if(!scene.Cast(refractedRay,refractionHit)){
						insideMaterial = false;
					}
				}else{

					refractedRay.origin = refractionHit.point + epsilon*currentNormal;
					insideMaterial = false;

				}

				

			}while(insideMaterial);

			refractedRay.generation = hit.ray.generation + 1;

			//now do refraction
			refractedColor = scene.Trace(refractedRay);

		}
	}


	//do the reflection
	//only do reflection if the reflectance is greater than zero
	//	this is an optimization so unnecessary reflections are not calculated
	if( (r.red + r.green + r.blue) > epsilon){
		reflectedColor = scene.Trace(reflectionRay);
	}
	
	//debug code
	//printf("Eye Vector: (%f,%f,%f)\n",E.x,E.y,E.z);
	//printf("Normal Vector: (%f,%f,%f)\n",N.x,N.y,N.z);
	//printf("Refraction Vector: (%f,%f,%f)\n",refractionDir.x,refractionDir.y,refractionDir.z);
	//printf("(N_1,N_2)=(%f,%f)\n\n",hit.ray.ref_index,refractedRay.ref_index);
	

	//now combine calculated color with reflected color
	finalColor = (-1*t + Color(1.0,1.0,1.0))*(colorWithLighting + r*reflectedColor) + t*refractedColor;

	return finalColor; 
    }

	/*calculate the refraction direction
	*	The derivation was inspired by Slides 22-23 of this lecture:
	*		http://graphics.ucsd.edu/courses/cse168_s06/ucsd/lecture03.pdf
	*/
Vec3 basic_shader::RefractionDirection(double n_1, double n_2, Vec3 incomingVector,Vec3 normalVector) const{

	double ratio = n_1/n_2; //ratio = 1.0/ratio;
	double cos_theta = normalVector*incomingVector;
	Vec3 Refrac_vertical = ratio*( normalVector*cos_theta - incomingVector);
	double cos_phi_squared = 1 - ( ratio*ratio * (1-cos_theta*cos_theta));
	
	//total internal reflection occurs
	if(cos_phi_squared < 0.0){
		return Vec3();
		//return Unit( ( 2.0 * ( incomingVector * normalVector ) ) * normalVector - incomingVector );
	}

	Vec3 Refrac_horizontal = -1*normalVector*sqrt(cos_phi_squared);
	Vec3 Refrac = Refrac_vertical + Refrac_horizontal;
	return Unit(Refrac);
}
