/***************************************************************************
* block.cpp    (primitive object plugin)                                   *
*                                                                          *
* The "Block" object defines an axis-aligned box.  Its constructor takes   *
* two vectors, which are taken to be the "min" and "max" corners of the    *
* Block (i.e. the three min coords, and the three max coords).             *
*                                                                          *
* History:                                                                 *
*   10/10/2004  Split off from objects.cpp file.                           *
*                                                                          *
***************************************************************************/
#include "toytracer.h"
#include "util.h"
#include "params.h"

struct Block : public Primitive {
    Block() {}
    Block( const Vec3 &Min, const Vec3 &Max );
    virtual bool Intersect( const Ray &ray, HitInfo &hitinfo ) const;
    virtual bool Inside( const Vec3 &P ) const;
    virtual Interval GetSlab( const Vec3 & ) const;
    virtual int GetSamples( const Vec3 &P, const Vec3 &N, Sample *samples, int n ) const;
    virtual Plugin *ReadString( const string &params );
    virtual string MyName() const { return "block"; }
    Vec3 Min; // Minimum coordinates along each axis.
    Vec3 Max; // Maximum coordinates along each axis.
    };

// Register the new object with the toytracer.  When this module is linked in, the 
// toytracer will automatically recognize the new objects and read them from sdf files.

REGISTER_PLUGIN( Block );

Block::Block( const Vec3 &Min_, const Vec3 &Max_ )
    {
    Min  = Min_;
    Max  = Max_;
    }

Plugin *Block::ReadString( const string &params ) // Read params from string.
    {
    ParamReader get( params );
    Vec3 Vmin, Vmax;
    if( get[MyName()] && get[Vmin] && get[Vmax] ) return new Block( Vmin, Vmax );
    return NULL;
    }

// Determine whether the given point is on or in the object.
bool Block::Inside( const Vec3 &P ) const
    {
    if( P.x < Min.x || P.x > Max.x ) return false;
    if( P.y < Min.y || P.y > Max.y ) return false;
    if( P.z < Min.z || P.z > Max.z ) return false;
    return true;
    }

Interval Block::GetSlab( const Vec3 &v ) const
    {
    const double
        ax = v.x * Min.x,
        ay = v.y * Min.y,
        az = v.z * Min.z,
        bx = v.x * Max.x,
        by = v.y * Max.y,
        bz = v.z * Max.z,
        a  = min( ax, bx ) + min( ay, by ) + min( az, bz ), 
        b  = max( ax, bx ) + max( ay, by ) + max( az, bz ); 
    return Interval( a, b ) / ( v * v );
    }

bool Block::Intersect( const Ray &ray, HitInfo &hitinfo ) const
    {
	/*
	* With this primitive, we need to compute the intersection
	*	with each axis.
	* This will let us get an interval of t values for x, y, and z
	*	The overlap of the intervals will be the regions of our box
	*	If the final interval is in our acceptable t range, the ray
	*		intersects the box
	* If the intervals can't overlap, then the ray does not intersect
	* 
	* Here is the procedure in more detail:
	*	(1) For x, compute t_at_xmin = (Min_x - origin_x)/(direction_x)
	*			   compute t_at_xmax = (Max_x - origin_x)/(direction_x)
	*	(2) if direction_x >= 0
	*				set tx_min = t_at_xmin; tx_max=t_at_xmax
	*		else
	*				set tx_max = t_at_xmin; tx_min=t_at_xmax
	*	(3) repeat (1) and (2) but for y and z
	*	(4) set (t_min,t_max) equal to the intersection of 
	*			(tx_min,tx_max), (ty_min,ty_max), (tz_min,tz_max)
	*/
	const double min_t = 1;
	const double max_t = 20;

	double t_at_xmin,t_at_xmax,t_at_ymin,t_at_ymax,t_at_zmin,t_at_zmax;
	double tmin, tmax, tymin, tymax, tzmin, tzmax;
	double dist;

	t_at_xmin = (Min.x - ray.origin.x)/(ray.direction.x); 
	t_at_xmax = (Max.x - ray.origin.x)/(ray.direction.x); 

	if(ray.direction.x >= 0){
		tmin = t_at_xmin;
		tmax = t_at_xmax;
	}else{
		tmin = t_at_xmax;
		tmax = t_at_xmin;
	}

	t_at_ymin = (Min.y - ray.origin.y)/(ray.direction.y); 
	t_at_ymax = (Max.y - ray.origin.y)/(ray.direction.y); 

	if(ray.direction.y >= 0){
		tymin = t_at_ymin;
		tymax = t_at_ymax;
	}else{
		tymin = t_at_ymax;
		tymax = t_at_ymin;
	}

	if( (tmin > tymax) || (tmax < tymin) ){
		return false;
	}

	if(tymin > tmin){
		tmin = tymin;
	}
	if(tymax < tmax){
		tmax = tymax;
	}

	t_at_zmin = (Min.z - ray.origin.z)/(ray.direction.z); 
	t_at_zmax = (Max.z - ray.origin.z)/(ray.direction.z); 

	if(ray.direction.z >= 0){
		tzmin = t_at_zmin;
		tzmax = t_at_zmax;
	}else{
		tzmin = t_at_zmax;
		tzmax = t_at_zmin;
	}

	if( (tmin > tzmax) || (tmax < tzmin) ){
		return false;
	}

	if(tzmin > tmin){
		tmin = tzmin;
	}
	if(tzmax < tmax){
		tmax = tzmax;
	}

	//look at tmin and tmax to see which one should be used
	if(tmin > 0){
		if(tmin > hitinfo.distance){
			return false;
		}
		dist = tmin;
	}else{
		//we are either inside the block or the block is behind us

		if(tmax <= 0){
			//the block is behind us
			return false;
		}
		if(tmax > hitinfo.distance){
			return false;
		}

		//we are definitely inside the block so set the intersection point to tmax
		dist = tmax;
	}
	
	/* Knowing dist will give us the face that we are looking at
	*	Since the faces are axis-aligned, the normals are
	*		<0,0,1> , <0,1,0>, <1,0,0> 
	*		<0,0,-1> , <0,-1,0>, or <-1,0,0>
	*		depending on which face it is
	*	We will see what face it is an use that for the normals
	*/
	Vec3 currentNormal = Vec3(0,0,0);
	if(dist == t_at_xmin){
		currentNormal = Vec3(-1,0,0);
	}else if(dist == t_at_xmax){
		currentNormal = Vec3(1,0,0);
	}else if(dist == t_at_ymin){
		currentNormal = Vec3(0,-1,0);
	}else if(dist == t_at_ymax){
		currentNormal = Vec3(0,1,0);
	}else if(dist == t_at_zmin){
		currentNormal = Vec3(0,0,-1);
	}else if(dist == t_at_zmax){
		currentNormal = Vec3(0,0,1);
	}

	hitinfo.distance = dist;
	hitinfo.point    = ray.origin + dist * ray.direction;
	hitinfo.object   = this;
	hitinfo.normal = currentNormal;
	return true;
    }

int Block::GetSamples( const Vec3 &P, const Vec3 &N, Sample *samples, int n ) const
    {
    // To be filled in later.
    return 0;
    }

