/***************************************************************************
* list.cpp   (aggregate object plugin)                                     *
*                                                                          *
* The "list" object is a trivial aggregate object.  It simply collects     *
* all children into a simple linked list and tests each ray against        *
* each object in the list.  It places a single bounding box around the     *
* entire collection, but not around each individual object.  This          *
* amounts to brute-force ray tracing.                                      *
*                                                                          *
* History:                                                                 *
*   10/16/2004  Changed the way the bounding box is computed.              *
*   10/16/2004  Added more documentation, and call to "Inverse" function.  *
*   10/06/2004  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#include "toytracer.h"
#include "util.h"
#include "params.h"

// Define the actual List object as a sub-class of the "Aggregate" class.
// This sub-class must define all the necessary virtual methods as well as
// any special data members that are specific to this type of object.
// (ALL access to this object will be through these virtual methods.)

struct List : public Aggregate { 
    List() { bbox = AABB::Null(); }
   ~List();
    virtual bool Intersect( const Ray &ray, HitInfo & ) const;
    virtual bool Inside( const Vec3 & ) const;
    virtual Interval GetSlab( const Vec3 & ) const;
    virtual Plugin *ReadString( const string &params );
    virtual string MyName() const { return "List"; }
    virtual void AddChild( Object * );
    virtual double Cost() const;
    AABB bbox;
    };

// Register the new object with the toytracer.  When this module is linked in, the 
// toytracer will automatically recognize the new objects and read them from sdf files.
REGISTER_PLUGIN( List );

List::~List()
    {
    // No need to do any memory deallocation here, as the destructor will
    // only be called when the image is complete (thus, we can let the system
    // do all the deallocation as a normal result of program termination).
    }

// Fill in all the virtual methods for the List object...

Plugin *List::ReadString( const string &params )
    {
    ParamReader get( params );
    // Every aggregate object should look for a "begin" preceeding its name.
    // An aggregate can also accept additional parameters on this line (e.g. the
    // size of a uniform grid, maximum depth of a BSP tree, etc.).  It is recommended
    // that missing parameters default to something reasonable.
    if( get["begin"] && get[MyName()] ) return new List();
    return NULL;
    }

bool List::Intersect( const Ray &ray, HitInfo &hitinfo ) const
    {
    // If the ray does not hit the bounding box, then the ray
    // misses the contents.
    if( !Hit( ray, bbox ) ) return false;

    bool found_a_hit = false;

    // Each intersector is ONLY allowed to write into the "HitInfo"
    // structure if it has determined that the ray hits the object
    // at a CLOSER distance than currently recorded in HitInfo.distance.
    // When a closer hit is found, the fields of the "HitInfo" structure
    // are updated to hold the material of the object that was just hit.
    // Return "true" if any object was hit.

    for( unsigned i = 0; i < NumChildren(); i++ )
        {
        const Object *obj = GetChild(i);
        if( obj != hitinfo.ignore && obj->Intersect( ray, hitinfo ) )
            {
            // We've hit an object, but continue looking as there might be a closer hit.
            found_a_hit = true; 
            }
        }
    // No need to fill in any fields of hitinfo, as the closest object
    // that was hit (if there was one) will have already done that.
    return found_a_hit;
    }

bool List::Inside( const Vec3 &P ) const
    {
    // If the point is not inside the bounding box of the list, then
    // there is no need to check the child objects.
    if( !::Inside( P, bbox ) ) return false;
    for( unsigned i = 0; i < NumChildren(); i++ )
        {
        const Object *obj = GetChild(i);
        if( obj->Inside( P ) ) return true;
        }
    return false;
    }

Interval List::GetSlab( const Vec3 &v ) const
    {
    Interval I = Interval::Null();
    for( unsigned i = 0; i < NumChildren(); i++ )
        {
        // Expand the interval for this slab to ensure that it
        // encloses the intervals enclosing each child object.
        I << GetChild(i)->GetSlab(v);
        }
    return I;
    }

double List::Cost() const
    {
    // Since the List object simply tests all its children, in the order they
    // were given, the cost of intersecting a ray with this type of object is
    // the sum of the costs of all the children.  (Hence, this is a VERY
    // inefficient aggregate object.)
    double cost = 0.0;
    for( unsigned i = 0; i < NumChildren(); i++ )
        cost += GetChild(i)->Cost();
    return cost;
    }

void List::AddChild( Object *obj )
    {
    // The list object supplies its own AddChild method so that it can
    // construct a bounding box containing all the child objects as they
    // are added.
    children.push_back( obj );
    bbox << GetBox( *obj );  // Grow the bounding box.
    }





