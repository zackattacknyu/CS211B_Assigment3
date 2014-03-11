/***************************************************************************
* torus.cpp   (primitive object plugin)                                    *
*                                                                          *
* The "torus" primitive object is situated on the x-y plane.  It has two   *
* parameters: the major radius (denoted by "a") and the minor radius       *
* (denoted by "b").  The ray-torus intersection test is performed          *
* analytically using a closed-form quartic polynomial root finder.         *
*                                                                          *
* History:                                                                 *
*   10/12/2005  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#include "toytracer.h"
#include "util.h"
#include "params.h"
#include "quartic.h"

struct torus : public Primitive {
    torus() {}
    torus( double major_radius, double minor_radius );
    virtual bool Intersect( const Ray &ray, HitInfo & ) const;
    virtual bool Inside( const Vec3 &P ) const; 
    virtual Interval GetSlab( const Vec3 & ) const;
    virtual int GetSamples( const Vec3 &P, const Vec3 &N, Sample *samples, int n ) const;
    virtual Plugin *ReadString( const string &params );
    virtual string MyName() const { return "torus"; }
    virtual double Cost() const { return 4.0; }
    double a;    // Major radius.
    double b;    // Minor radius.
    double a2b2; // a * a + b * b
    double rad;  // a + b
    double rad2; // rad * rad 
    AABB   bbox;
    };

REGISTER_PLUGIN( torus );

torus::torus( double major_radius, double minor_radius )
    {
    // Pre-compute the often-used values that depend only on the parameters.
    a    = major_radius;
    b    = minor_radius;
    a2b2 = a * a + b * b;
    rad  = a + b;
    rad2 = rad * rad;
    bbox = AABB( Interval( -rad, rad ), Interval( -rad, rad ), Interval( -b, b ) );
    }

Plugin *torus::ReadString( const string &params ) // Reads params from a string.
    {
    double ra, rb;
    ParamReader get( params );
    if( get[MyName()] && get[ra] && get[rb] ) return new torus( ra, rb );
    return NULL;
    }

Interval torus::GetSlab( const Vec3 &v ) const
    {
    double d = a * ( v.x * v.x + v.y * v.y ) / ( v * v ) + b / Length(v);
    return Interval( -d, d );
    }

bool torus::Inside( const Vec3 &P ) const
    {
    // Quick up-down reject test.
    if( P.z > b || P.z < -b ) return false;
    const double r = sqrt( P.x * P.x + P.y * P.y );
    // Quick radial reject test.
    if( r > rad || r < a - b ) return false;
    const double scale = a / r;
    Vec3 Q( scale * P.x, scale * P.y, 0.0 );
    return dist( P, Q ) <= b;
    }

bool torus::Intersect( const Ray &ray, HitInfo &hitinfo ) const
    {
    // THINGS TO DO: Should reposition the ray origin to the intersection
    // with the bounding box.  This will increase the accuracy of the root
    // finder.
    roots rts;
    const Vec3 R( ray.direction );
    const Vec3 Q( ray.origin    );

    // First look for easy reject cases.
    // Test the bounding planes on the top and bottom: z = b and z = -b.
    if( R.z >= 0.0 )
        {
        // Above the torus looking up.
        if( Q.z >=  b ) return false;
        }
    else
        {
        // Below the torus looking down.
        if( Q.z <= -b ) return false;
        }

    // See if the ray trivially misses the infinite bounding cylinder.
    // If not, the results of this computation will be used for root finding.
    const double qq = Q.x * Q.x + Q.y * Q.y;
    const double qr = Q.x * R.x + Q.y * R.y;
    if( qr > 0.0 && qq >= rad2 ) return false;

    // Perform a bounding-box test.
    if( !Hit( ray, bbox ) ) return false;

    const double QQ = qq + Q.z * Q.z;
    const double QR = qr + Q.z * R.z;
    const double a2 = a * a;

    // Set up and solve the quartic polynomial that results from plugging
    // the ray equation into the implicit equation of the torus.  That is,
    //
    //   || P  -  a Q(P) || = b    (implicit euqation of torus.)
    //
    //    P  =  Q + R s            (parametric form of a ray.)
    //
    // where Q(P) is the closest point to P that is on the unit-radius circle,
    // centered at the origin, in the x-y plane.  The function Q is given by
    //
    //    Q(P)  =  ( Px, Py, 0 ) / sqrt( Px Px  +  Py Py )
    //
    // When the radical are removed via squaring, the result is a quartic
    // equation in s, the distance along the ray to the point of intersection.
    //
    SolveQuartic(
        sqr( QQ - a2b2 ) - 4.0 * a2 * ( b * b - Q.z * Q.z ),      // Constant
        4.0 * QR * ( QQ - a2b2 ) + 8.0 * a2 * Q.z * R.z,          // Linear
        2.0 * ( QQ - a2b2 ) + 4.0 * ( QR * QR + a2 * R.z * R.z ), // Quadratic
        4.0 * QR,                                                 // Cubic
        1.0,                                                      // Quartic
        rts // The real roots returned by the solver.
        );

    // If the smallest positive root is larger than the distance already
    // in hitinfo, then we do not consider it a hit.
    double s = hitinfo.distance;
    if( !rts.MinPositiveRoot( s ) ) return false;

    // Compute the actual point of intersection using the distance.
    Vec3 P( ray.origin + s * R );

    hitinfo.distance = s;
    hitinfo.point    = P;
    hitinfo.normal   = Unit( P - a * Unit( P.x, P.y, 0.0 ) ); 
    hitinfo.object   = this;
    return true;
    }

int torus::GetSamples( const Vec3 &P, const Vec3 &N, Sample *samples, int n ) const
    {
    // To be supplied...
    return 0;
    }

