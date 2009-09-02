#include "scalarfield.hpp"
#include <iostream>
#include <string.h>
#include <cmath>
#include <limits>


ScalarField::ScalarField( const Geometry &g )
    : _geom_mode(g.geom_mode()), _size(g.size()), _origo(g.origo()), _h(g.h())
{
    // Calculate vector max
    _max = Vec3D( _origo(0)+_h*(_size[0]-1),
		  _origo(1)+_h*(_size[1]-1),
		  _origo(2)+_h*(_size[2]-1) );
    _div_h = 1.0/_h;
    _F = new double[_size[0]*_size[1]*_size[2]];
    memset( _F, 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
}


ScalarField::ScalarField( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h )
    : _geom_mode(geom_mode), _size(size), _origo(origo), _h(h)
{
    if( _geom_mode == MODE_3D ) {
	if( _size[0] < 3 || _size[1] < 3 || _size[2] < 3 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    } else if( _geom_mode == MODE_2D || _geom_mode == MODE_CYL ) {
	if( _size[0] < 3 || _size[1] < 3 || _size[2] != 1 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    } else {
	if( _size[0] < 3 || _size[1] != 1 || _size[2] != 1 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    }

    // Calculate vector max
    _max = Vec3D( _origo(0)+_h*(_size[0]-1),
		  _origo(1)+_h*(_size[1]-1),
		  _origo(2)+_h*(_size[2]-1) );
    _div_h = 1.0/_h;
    _F = new double[_size[0]*_size[1]*_size[2]];
    memset( _F, 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
}


ScalarField::ScalarField( std::istream &s )
{
    _geom_mode = (geom_mode_e)read_int32( s );
    _size      = Int3D( s );
    _origo     = Vec3D( s );
    _h         = read_double( s );
    // Calculate vector max
    _max = Vec3D( _origo(0)+_h*(_size[0]-1),
		  _origo(1)+_h*(_size[1]-1),
		  _origo(2)+_h*(_size[2]-1) );
    _div_h     = 1.0/_h;
    _F = new double[_size[0]*_size[1]*_size[2]];
    read_compressed_block( s, _size[0]*_size[1]*_size[2]*sizeof(double), (int8_t *)_F );
}


ScalarField::ScalarField( const ScalarField &f )
    : _geom_mode(f._geom_mode), _size(f._size), _origo(f._origo), _max(f._max), _h(f._h), 
      _div_h(f._div_h)
{
    _F = new double[_size[0]*_size[1]*_size[2]];
    memcpy( _F, f._F, _size[0]*_size[1]*_size[2]*sizeof(double) );
}


ScalarField::~ScalarField()
{
    delete [] _F;
}


void ScalarField::clear()
{
    memset( _F, 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
}


void ScalarField::reset( geom_mode_e geom_mode, Int3D size, 
			 Vec3D origo, double h )
{
    _geom_mode = geom_mode;
    _size      = size;
    _origo     = origo;
    _h         = h;
    _max = Vec3D( _origo(0)+_h*(_size[0]-1),
		  _origo(1)+_h*(_size[1]-1),
		  _origo(2)+_h*(_size[2]-1) );
    _div_h     = 1.0/h;
    delete [] _F;
    _F = new double[_size[0]*_size[1]*_size[2]];
    memset( _F, 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
}


void ScalarField::get_minmax( double &min, double &max ) const
{
    min = _F[0];
    max = _F[0];

    size_t ncount = _size[0]*_size[1]*_size[2];
    std::cout << "ncount = " << ncount << "\n";
    for( size_t a = 1; a < ncount; a++ ) {
	if( _F[a] < min )
	    min = _F[a];
	if( _F[a] > max )
	    max = _F[a];
	std::cout << "a = " << a << ": " 
		  << "min = " << min << " "
		  << "max = " << max << "\n";
    }
}


void ScalarField::epot_get_minmax( const Geometry &g, double &min, double &max ) const
{
    min = std::numeric_limits<double>::infinity();
    max = -std::numeric_limits<double>::infinity();

    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t a = 0; a < ncount; a++ ) {
	if( g.mesh( a ) <= -7 )
	    // Skip electrode edges
	    continue; 
	if( _F[a] < min )
	    min = _F[a];
	if( _F[a] > max )
	    max = _F[a];
    }
}


ScalarField &ScalarField::operator=( const ScalarField &f )
{
    _geom_mode = f._geom_mode;
    _size      = f._size;
    _origo     = f._origo;
    _h         = f._h;
    _div_h     = f._div_h;
    delete [] _F;
    _F = new double[_size[0]*_size[1]*_size[2]];
    memcpy( _F, f._F, _size[0]*_size[1]*_size[2]*sizeof(double) );
    return( *this );
}


ScalarField &ScalarField::operator+=( const ScalarField &f )
{
    if( _geom_mode != f._geom_mode || _size != f._size ||
	_origo != f._origo || _h != f._h )
	throw( Error( ERROR_LOCATION, "non-matching fields" ) );
    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t a = 0; a < ncount; a++ )
	_F[a] += f._F[a];
    return( *this );
}


ScalarField &ScalarField::operator*=( double x )
{
    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t a = 0; a < ncount; a++ )
	_F[a] *= x;
    return( *this );
}


ScalarField &ScalarField::operator/=( double x )
{
    double xi = 1.0/x;
    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t a = 0; a < ncount; a++ )
	_F[a] *= xi;
    return( *this );
}


double ScalarField::operator()( Vec3D x ) const
{
    switch( _geom_mode ) {
    case MODE_1D:
    {
	int32_t i = (int32_t)floor( (x[0]-_origo[0])*_div_h );
	if( i < 0 )
	    i = 0;
	else if( i >= _size[0]-1 )
	    i = _size[0]-2;

	double t = _div_h*( x[0]-(i*_h+_origo[0]) );

	return( (1.0-t)*_F[i] + t*_F[i+1] );
	break;
    }
    case MODE_2D:
    case MODE_CYL:
    {
	int32_t i = (int32_t)floor( (x[0]-_origo[0])*_div_h );
	int32_t j = (int32_t)floor( (x[1]-_origo[1])*_div_h );
	if( i < 0 )
	    i = 0;
	else if( i >= _size[0]-1 )
	    i = _size[0]-2;
	if( j < 0 )
	    j = 0;
	else if( j >= _size[1]-1 )
	    j = _size[1]-2;

	double t = _div_h*( x[0]-(i*_h+_origo[0]) );
	double u = _div_h*( x[1]-(j*_h+_origo[1]) );

	int32_t ptr = _size[0]*j + i;
	return( (1.0-t)*(1.0-u)*_F[ptr] +
		(    t)*(1.0-u)*_F[ptr+1] +
		(1.0-t)*(    u)*_F[ptr+_size[0]] +
		(    t)*(    u)*_F[ptr+1+_size[0]] );
	break;
    }
    default:
    {
	int32_t i = (int32_t)floor( (x[0]-_origo[0])*_div_h );
	int32_t j = (int32_t)floor( (x[1]-_origo[1])*_div_h );
	int32_t k = (int32_t)floor( (x[2]-_origo[2])*_div_h );
	if( i < 0 )
	    i = 0;
	else if( i >= _size[0]-1 )
	    i = _size[0]-2;
	if( j < 0 )
	    j = 0;
	else if( j >= _size[1]-1 )
	    j = _size[1]-2;
	if( k < 0 )
	    k = 0;
	else if( k >= _size[2]-1 )
	    k = _size[2]-2;

	double t = _div_h*( x[0]-(i*_h+_origo[0]) );
	double u = _div_h*( x[1]-(j*_h+_origo[1]) );
	double v = _div_h*( x[2]-(k*_h+_origo[2]) );

	int32_t b  = _size[0]*_size[1];
	int32_t ptr = b*k + _size[0]*j + i;
	return( (1.0-t)*(1.0-u)*(1.0-v)*_F[ptr] +
		(    t)*(1.0-u)*(1.0-v)*_F[ptr+1] +
		(1.0-t)*(    u)*(1.0-v)*_F[ptr+_size[0]] +
		(    t)*(    u)*(1.0-v)*_F[ptr+1+_size[0]] +
		(1.0-t)*(1.0-u)*(    v)*_F[ptr+b] +
		(    t)*(1.0-u)*(    v)*_F[ptr+1+b] +
		(1.0-t)*(    u)*(    v)*_F[ptr+_size[0]+b] +
		(    t)*(    u)*(    v)*_F[ptr+1+_size[0]+b] );
	break;
    }
    }

}


void ScalarField::save( std::ostream &s ) const
{
    write_int32( s, _geom_mode );
    _size.save( s );
    _origo.save( s );
    write_double( s, _h );
    write_compressed_block( s, _size[0]*_size[1]*_size[2]*sizeof(double), (int8_t *)_F );
}


void ScalarField::debug_print( void ) const
{
    int32_t a;
    std::cout << "**ScalarField\n";
    if( _geom_mode == MODE_1D )
	std::cout << "geom_mode = MODE_1D\n";
    else if( _geom_mode == MODE_2D )
	std::cout << "geom_mode = MODE_2D\n";
    else if( _geom_mode == MODE_CYL )
	std::cout << "geom_mode = MODE_CYL\n";
    else if( _geom_mode == MODE_3D )
	std::cout << "geom_mode = MODE_3D\n";
    else
	std::cout << "geom_mode = Unknown\n";
    std::cout << "size = (" 
	      << _size[0] << ", "
	      << _size[1] << ", "
	      << _size[2] << ")\n";
    std::cout << "origo = (" 
	      << _origo[0] << ", "
	      << _origo[1] << ", "
	      << _origo[2] << ")\n";
    std::cout << "h = " << _h << "\n";
    std::cout << "div_h = " << _div_h << "\n";
    std::cout << "F = (";
    for( a = 0; a < _size[0]*_size[1]*_size[2]-1; a++ )
	std::cout << _F[a] << ", ";
    if( a < _size[0]*_size[1]*_size[2] )
	std::cout << _F[a] << ")\n";

}















