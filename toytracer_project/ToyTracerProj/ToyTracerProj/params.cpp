/***************************************************************************
* params.cpp                                                               *
*                                                                          *
* A tool for reading parameter strings from sdf files.  Each method        *
* looks for some specific pattern at the beginning of a string and does    *
* three things if it is found:                                             *
*    1) The pattern is removed from the string.                            *
*    2) The pattern is parsed and assigned to the argument, if appropriate.*
*    3) True is returned as the function value.                            *
* If the pattern is not found, the string is left intact and False is      *
* returned as the function value.                                          *
*                                                                          *
* History:                                                                 *
*   10/04/2005  Initial coding.                                            *
*                                                                          *
***************************************************************************/
#include "params.h"

void ParamReader::SkipBlanks()
    {
    while( params[0] == ' ' || params[0] == '\t' ) params.erase(0,1);
    }

bool ParamReader::operator[]( string field )
    {
    SkipBlanks();
    if( params.find( field, 0 ) == 0 )
        {
        params.erase( 0, field.length() );
        return true;
        }
    return false;
    }

bool ParamReader::operator[]( Color &c )
    {
    if( sscanf( params.c_str(), " [ %lf , %lf , %lf ]", &c.red, &c.green, &c.blue ) == 3 )
        {
        int len = 1 + params.find( "]", 0 );
        params.erase( 0, len );
        return true;
        }
    return false;
    }

bool ParamReader::operator[]( Vec3 &v )
    {
    if( sscanf( params.c_str(), " ( %lf , %lf , %lf )", &v.x, &v.y, &v.z ) == 3 )
        {
        int len = 1 + params.find( ")", 0 );
        params.erase( 0, len );
        return true;
        }
    return false;
    }

bool ParamReader::operator[]( Vec2 &v )
    {
    if( sscanf( params.c_str(), " ( %lf , %lf )", &v.x, &v.y ) == 2 )
        {
        int len = 1 + params.find( ")", 0 );
        params.erase( 0, len );
        return true;
        }
    return false;
    }

bool ParamReader::operator[]( Interval &I )
    {
    if( sscanf( params.c_str(), " ( %lf , %lf )", &I.min, &I.max ) == 2 )
        {
        int len = 1 + params.find( ")", 0 );
        params.erase( 0, len );
        return true;
        }
    return false;
    }

bool ParamReader::operator[]( double &x )
    {
    SkipBlanks();
    if( sscanf( params.c_str(), "%lf", &x ) == 1 )
        {
        int len = 1 + params.find( " ", 0 );
        if( len > 1 ) params.erase( 0, len );
        return true;
        }
    return false;
    }

bool ParamReader::operator[]( unsigned &x )
    {
    SkipBlanks();
    if( sscanf( params.c_str(), "%u", &x ) == 1 )
        {
        int len = 1 + params.find( " ", 0 );
        if( len > 1 ) params.erase( 0, len );
        return true;
        }
    return false;
    }

bool ParamReader::operator[]( Mat3x4 &M )
    {
    // Read the given string followed by a 3x4 matrix, in row order, enclosed in parens.
    int count = sscanf(
        params.c_str(),
        " ( %lf , %lf , %lf , %lf ; %lf , %lf , %lf , %lf ; %lf , %lf , %lf , %lf ) ",
        &M.mat(0,0), &M.mat(0,1), &M.mat(0,2), &M.vec.x, 
        &M.mat(1,0), &M.mat(1,1), &M.mat(1,2), &M.vec.y,
        &M.mat(2,0), &M.mat(2,1), &M.mat(2,2), &M.vec.z
        );
    if( count == 12 )
        {
        int len = 1 + params.find( ")", 0 );
        params.erase( 0, len );
        return true;
        }
    return false;
    }

