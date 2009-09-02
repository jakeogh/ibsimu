#include <iostream>
#include "csgobject_solid.hpp"
#include "error.hpp"


CSGObjectSolid::CSGObjectSolid( std::istream &s )
{
    _object = 0;
}


bool CSGObjectSolid::inside( const Vec3D &x ) const
{
    return( _object->inside( Vector4( x[0], x[1], x[2], 1.0 ) ) );
}


void CSGObjectSolid::debug_print( void ) const
{

}


void CSGObjectSolid::save( std::ostream &fout ) const
{
    write_int32( fout, FILEID_CSGSOLID );
}















