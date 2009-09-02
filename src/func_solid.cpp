#include <iostream>
#include "func_solid.hpp"
#include "error.hpp"


bool FuncSolid::inside( const Vec3D &x ) const
{
    if( !_func )
	throw( Error( ERROR_LOCATION, "solid function not defined" ) );
    return( _func( x[0], x[1], x[2] ) );
}


FuncSolid::FuncSolid( std::istream &s )
{
    _func = 0;
}


void FuncSolid::debug_print( void ) const
{
    std::cout << "**FuncSolid\n";
    std::cout << "func = " << _func << "\n";
}


void FuncSolid::save( std::ostream &fout ) const
{
    write_int32( fout, FILEID_FUNCSOLID );
}













