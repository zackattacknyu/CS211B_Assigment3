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
	virtual Color basic_shader::GetDiffuseColor(Vec3 lightVector,Vec3 normal,Color emission, Color diffuse) const;
    };

static const double attenuation_a = 0.0;
static const double attenuation_b = 0.0;
static const double attenuation_c = 0.02;

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
	double shadowFactor = 0;
	bool objectWasHit = false;

	int numGridVals = 30;
	float numSamples = 5;
	int totalNumGridVals = numGridVals*numGridVals*numSamples;
	int totalArrayValues = totalNumGridVals*2;
	float totalGridVals = numGridVals;
	float currentXvalue;
	float currentYvalue;
	float interval = 2.0f/totalGridVals;
	float currentDist;
	bool posZinHemisphere,negZinHemisphere;
	int currentInd1,currentInd2;
	float randomX,randomY;
	float currentPosZvalue,currentNegZvalue;
	//Vec3 currentNormal = Vec3(0.0f,0.0f,1.0f);
	Vec3 currentNormal = N;
	float dotProdPosZ,dotProdNegZ,dotProdXYpart;
	Vec3 positiveZvector;
	Vec3 negativeZvector;
	int numLightsHit = 0;
	
	/*
	Surface Area of sphere from each patch P is approximately area(P)*(1/z') where z' is taken at the center
	This comes from the fact that the surface area is integral_P (1/z)
	*/
	float minZvalue=sqrt(interval)*sqrt(2-interval);;
	
	//used to calculate surface area
	float centerXvalue,centerYvalue,centerZvalue,approxSurfaceArea;

	//temp variable. default emission of the light blocks
	Color defaultEmission = Color(1.0,1.0,1.0);
	double emissionFactor = 30.0;

	Color posZpatchValue = Color();
	Color negZpatchValue = Color();
	Color posZpatchSpecValue = Color();
	Color negZpatchSpecValue = Color();
	double radius,patchFormFactor;
	Vec3 otherNormal;
	float currentEnergy = 0;

	for(int xInd = 0; xInd < numGridVals; xInd++){
		for(int yInd = 0; yInd < numGridVals; yInd++){

			centerXvalue = -1 + interval*xInd + interval*0.5;
			centerYvalue = -1 + interval*yInd + interval*0.5;
			currentDist = centerXvalue*centerXvalue + centerYvalue*centerYvalue;
			centerZvalue = sqrt(1-currentDist);

			
			if(centerZvalue < minZvalue){
				centerZvalue = minZvalue;
			}

			approxSurfaceArea = interval*interval*(1.0f/centerZvalue);	
			

			posZpatchValue = Color();
			negZpatchValue = Color();

			for(int sampleInd = 0; sampleInd < numSamples; sampleInd++){

				randomX = (double)rand() / RAND_MAX;
				randomY = (double)rand() / RAND_MAX;

				currentXvalue = -1+interval*xInd + interval*randomX;
				currentYvalue = -1+interval*yInd + interval*randomY;
				
				currentDist = currentXvalue*currentXvalue + currentYvalue*currentYvalue;

				posZinHemisphere = false;
				negZinHemisphere = false;

				currentInd1 = (xInd*numGridVals + yInd)*numSamples + sampleInd;
				currentInd2 = currentInd1 + totalNumGridVals;

				//makes sure it is inside the unit disk
				if(currentDist <= 1){

					currentPosZvalue = sqrt(1-currentDist);
					currentNegZvalue = -currentPosZvalue;

					dotProdXYpart = currentNormal.x*currentXvalue + currentNormal.y*currentYvalue;
					dotProdPosZ = currentNormal.z*currentPosZvalue + dotProdXYpart;
					dotProdNegZ = currentNormal.z*currentNegZvalue + dotProdXYpart;

					posZinHemisphere = (dotProdPosZ >=0);
					negZinHemisphere = (dotProdNegZ >=0);

				}


				if(posZinHemisphere){

					positiveZvector = Vec3(currentXvalue,currentYvalue,currentPosZvalue);

					//light ray to case to determine occulsion
					ray.origin = P;
					ray.direction = positiveZvector;
					HitInfo objectHit;
					objectHit.distance = Infinity;

					shadowFactor = 1;
					if(scene.Cast(ray,objectHit) ){
						if(objectHit.object != NULL){
							Vec3 LightPos = objectHit.point;
							
							if(LightPos.z > 12.0){
								lightVector = LightPos - P;
								numLightsHit = numLightsHit + 1;
								radius = Length(lightVector);
								otherNormal = Unit(objectHit.normal);
								lightVector = Unit(lightVector);

								//this is the current patches approximation of F_ij
								//patchFormFactor = ((abs(lightVector*currentNormal))*(abs(-1*lightVector*otherNormal)))/(Pi*radius*radius);
								patchFormFactor = ((lightVector*currentNormal)*(-1*lightVector*otherNormal))/(Pi*radius*radius);
								//printf("radius: %f\n",radius);
								currentEnergy = currentEnergy + emissionFactor*patchFormFactor;
								//posZpatchValue = posZpatchValue + GetDiffuseColor(lightVector,N,defaultEmission,diffuse);
							}
						}
					}

				}

				if(negZinHemisphere){
					negativeZvector = Vec3(currentXvalue,currentYvalue,currentNegZvalue);
					//printf("(x,y,z)=(%f,%f,%f)\n",gridXvalues[currentInd2],gridYvalues[currentInd2],gridZvalues[currentInd2]);

					//light ray to case to determine occulsion
					ray.origin = P;
					ray.direction = negativeZvector;
					HitInfo objectHit;
					objectHit.distance = Infinity;

					shadowFactor = 1;
					if(scene.Cast(ray,objectHit) ){
						if(objectHit.object != NULL){
							Vec3 LightPos = objectHit.point;
							
							if(LightPos.z > 12.0){
								lightVector = LightPos - P;
								numLightsHit = numLightsHit + 1;
								radius = Length(lightVector);
								
								otherNormal = Unit(objectHit.normal);
								lightVector = Unit(lightVector);

								//this is the current patches approximation of F_ij
								//patchFormFactor = ((abs(lightVector*currentNormal))*(abs(-1*lightVector*otherNormal)))/(Pi*radius*radius);
								patchFormFactor = ((lightVector*currentNormal)*(-1*lightVector*otherNormal))/(Pi*radius*radius);
								currentEnergy = currentEnergy + emissionFactor*patchFormFactor;
								
								//negZpatchValue = negZpatchValue + GetDiffuseColor(lightVector,N,defaultEmission,diffuse);
							}
						}
					}
				}

			}

			//diffuseColor = diffuseColor + (posZpatchValue/numSamples)*approxSurfaceArea;
			//diffuseColor = diffuseColor + (negZpatchValue/numSamples)*approxSurfaceArea;
			//printf("Approx Surface Area:%f\n",approxSurfaceArea);

			//diffuseColor = diffuseColor + posZpatchValue;
			//diffuseColor = diffuseColor + negZpatchValue;

			if(diffuseColor.blue > 0.5 || diffuseColor.green > 0.5 || diffuseColor.red > 0.5){
				//printf("Current Diffuse Color: (%f,%f,%f)\n",diffuseColor.blue,diffuseColor.green,diffuseColor.red);
			}
			
			
		}
	}
	
	//makes sure to include the light vectors in the calculations
	for( unsigned i = 0; i < scene.NumLights(); i++ )
        {
        const Object *light = scene.GetLight(i);
        Color emission = light->material->emission;
        AABB box = GetBox( *light );
        Vec3 LightPos( Center( box ) ); 

		//gets the light Vector
		lightVector = LightPos - P;
		lightDistance = Length(lightVector);
		Vec3 unitLightVector = Unit(lightVector);
		
		//gets the attenuation factor
		attenuation = 1/(attenuation_a + attenuation_b*lightDistance + attenuation_c*lightDistance*lightDistance);

		float dotProd = currentNormal.x*unitLightVector.x + currentNormal.y*unitLightVector.y + currentNormal.z*unitLightVector.z;

		//light vector in unit hemipshere
		if(dotProd >= 0){

			//light ray to case to determine occulsion
			ray.origin = P;
			ray.direction = unitLightVector;
			HitInfo objectHit;
			objectHit.distance = Infinity;

			if(scene.Cast(ray,objectHit) ){
				if(objectHit.object != NULL){
					Vec3 LightPos = objectHit.point;
							
					if(LightPos.z > 12.0){
						numLightsHit = numLightsHit + 1;
						lightVector = LightPos - P;
						//diffuseColor = diffuseColor + GetDiffuseColor(lightVector,N,defaultEmission,diffuse);

						radius = Length(lightVector);
								
						otherNormal = Unit(objectHit.normal);
						lightVector = Unit(lightVector);

						//this is the current patches approximation of F_ij
						//patchFormFactor = ((abs(lightVector*currentNormal))*(abs(-1*lightVector*otherNormal)))/(Pi*radius*radius);
						patchFormFactor = ((lightVector*currentNormal)*(-1*lightVector*otherNormal))/(Pi*radius*radius);
						currentEnergy = currentEnergy + emissionFactor*patchFormFactor;
					}
				}
			}

		}

		

		
    }

	diffuseColor = currentEnergy*defaultEmission;

	

	float numLights = numLightsHit;
	//diffuseColor = diffuseColor/(numLights*Pi);
	if(numLightsHit > 1){
		//printf("Number of Lights Hit:%d\n",numLightsHit);
		//printf("Original Color: (%f,%f,%f)\n",diffuse.blue,diffuse.green,diffuse.red);
		//printf("Diffuse Color: (%f,%f,%f)\n\n",diffuseColor.blue,diffuseColor.green,diffuseColor.red);

	}

	//colorWithLighting = color + diffuseColor*diffuse + specularColor*specular;
	//colorWithLighting = diffuseColor*diffuse;
	colorWithLighting = diffuseColor*diffuse + diffuse*0.4;

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
	
	//printf("diffuseColor: (%f,%f,%f)\n",colorWithLighting.red,colorWithLighting.green,colorWithLighting.blue);

	//now combine calculated color with reflected color
	//finalColor = (-1*t + Color(1.0,1.0,1.0))*(colorWithLighting + r*reflectedColor) + t*refractedColor;

	return colorWithLighting;
	//return finalColor; 
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

Color basic_shader::GetDiffuseColor(Vec3 lightVector,Vec3 normal,Color emission, Color diffuse) const{

	
	//printf("Hit Point: (%f,%f,%f)\n",LightPos.x,LightPos.y,LightPos.z);

	//gets the light Vector
	float lightDistance = Length(lightVector);

	lightVector = Unit(lightVector);

	//gets the attenuation factor
	//double attenuation = 1/(attenuation_a + attenuation_b*lightDistance + attenuation_c*lightDistance*lightDistance);
	double attenuation = 1;
	//gets the diffuse component
	float diffuseFactor = max(0,lightVector*normal);
		
	//gets the specular component
	//currentR = Unit(2.0*(N*lightVector)*N - lightVector);
	//specularFactor = max(0, pow(currentR*E,e) );

	if(diffuseFactor == 0){
		//specularFactor = 0;
	}

	//calculate the new color
	//specularColor = specularColor + shadowFactor*(attenuation*specularFactor)*emission;
	return (attenuation*diffuseFactor)*emission;

}
