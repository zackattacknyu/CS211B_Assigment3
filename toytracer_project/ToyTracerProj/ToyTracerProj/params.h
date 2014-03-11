/***************************************************************************
* params.h                                                                 *
*                                                                          *
* A tool for reading parameter strings from sdf files.                     *
*                                                                          *
* History:                                                                 *
*   10/04/2005  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#ifndef __PARAMS_INCLUDED__
#define __PARAMS_INCLUDED__

#include "toytracer.h"

// The ParamReader class makes the ReadString method of each Plugin very
// easy to write.  Given the parameters as a string, each "[]" operator
// attempts to strip off the given item, and returns "true" if successful.

class ParamReader {
    public:
        ParamReader( string s ) { params = s; }
        bool operator[]( string field );
        bool operator[]( Color    & );
        bool operator[]( Vec3     & );
        bool operator[]( Vec2     & );
        bool operator[]( double   & );
        bool operator[]( unsigned & );
        bool operator[]( Interval & ); 
        bool operator[]( Mat3x4   & );
    private:
        void SkipBlanks();
        string params;
	};

#endif

