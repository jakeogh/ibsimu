/*! \file vectorfield.cpp
 *  \brief Source code for vectorfield.cpp
 */

/* Copyright (c) 2005-2010 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * taneli.kalvas@jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */

#include <limits>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include "vectorfield.hpp"
#include "ibsimu.hpp"


VectorField::VectorField( const Geometry &g, bool fout[3] )
    : _geom_mode(g.geom_mode()), _size(g.size()), _origo(g.origo()), _h(g.h())
{
    // Check mesh size legality
    if( _size[0] < 1 || _size[1] < 1 || _size[2] < 1 )
	throw( Error( ERROR_LOCATION, "illegal mesh size" ) );

    _div_h = 1.0/_h;
    for( size_t i = 0; i < 3; i++ ) {
	if( fout[i] ) {
	    _F[i] = new double[_size[0]*_size[1]*_size[2]];
	    memset( _F[i], 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
	} else {
	    _F[i] = NULL;
	}
    }
}


VectorField::VectorField( geom_mode_e geom_mode, bool fout[3], Int3D size, 
			  Vec3D origo, double h )
    : _geom_mode(geom_mode), _size(size), _origo(origo), _h(h)
{
    // Check mesh size legality
    if( _size[0] < 1 || _size[1] < 1 || _size[2] < 1 )
	throw( Error( ERROR_LOCATION, "illegal mesh size" ) );

    _div_h = 1.0/_h;
    for( size_t i = 0; i < 3; i++ ) {
	if( fout[i] ) {
	    _F[i] = new double[_size[0]*_size[1]*_size[2]];
	    memset( _F[i], 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
	} else {
	    _F[i] = NULL;
	}
    }
}


bool VectorField::parse_line( const std::string &str, double c[6], double xscale, double fscale, 
			      size_t cdim, size_t fdim, const std::string &filename, size_t linec )
{
    const char *ptr = str.c_str();

    // Skip leading white space
    while( isspace(*ptr) ) ptr++;

    // Check if line contained only white space
    if( *ptr == '\n' || *ptr == '\r' || *ptr == '\0' )
	return( false );

    // Read cdim coordinates
    for( size_t i = 0; i < cdim; i++ ) {
	while( isspace(*ptr) ) ptr++;
	if( !isdigit(*ptr) && *ptr != '-' && *ptr != '+' && *ptr != '.' ) {
	    throw( Error( ERROR_LOCATION, "unexpected input reading coordinate data in " + 
			  filename + ":" + to_string( linec ) ) );
	}
	c[i] = xscale*strtod( ptr, (char **)(&ptr) );
    }

    // Read fdim field values
    for( size_t i = 0; i < fdim; i++ ) {
	while( isspace(*ptr) ) ptr++;
	if( !isdigit(*ptr) && *ptr != '-' && *ptr != '+' && *ptr != '.' ) {
	    throw( Error( ERROR_LOCATION, "unexpected input reading field data in " +
			  filename + ":" + to_string( linec ) ) );
	}
	c[3+i] = fscale*strtod( ptr, (char **)(&ptr) );
    }

    // Skip tailing white space
    while( isspace(*ptr) ) ptr++;

    // Check if line contains extra input
    if( *ptr != '\n' && *ptr != '\r' && *ptr != '\0' )
	throw( Error( ERROR_LOCATION, "unexpected input after reading field data in " +
		      filename + ":" + to_string( linec ) ) );

    return( true );
}


VectorField::VectorField( geom_mode_e geom_mode, bool fout[3], double xscale, 
			  double fscale, std::string filename )
{
    if( ibsimu.get_verbose_output() )
	std::cout << "Reading vector field from " << filename << "\n";
    
    // Read through data to find out mesh parameters
    size_t i;
    Int3D size;
    Vec3D origo( std::numeric_limits<double>::infinity(), 
		 std::numeric_limits<double>::infinity(), 
		 std::numeric_limits<double>::infinity() );
    Vec3D max( -std::numeric_limits<double>::infinity(), 
	       -std::numeric_limits<double>::infinity(), 
	       -std::numeric_limits<double>::infinity() );
    double h = 0.0;

    // Set number of dimensions (cdim) 
    size_t cdim;
    switch( geom_mode ) {
    case MODE_3D:
	cdim = 3;
	break;
    case MODE_2D:
	cdim = 2;
	break;
    case MODE_CYL:
	cdim = 2;
	break;
    case MODE_1D:
	cdim = 1;
	break;
    default:
	throw( Error( ERROR_LOCATION, "unsupported geometry mode" ) );
	break;
    }

    // Set number of field components (fdim)
    int fdim = 0;
    for( int a = 0; a < 3; a++ ) {
	if( fout[a] )
	    fdim++;
    }
    if( fdim == 0 )
	throw( Error( ERROR_LOCATION, "no field components to read" ) );

    double c[6], c2[6];
    int rec = 0;

    std::ifstream fin( filename.c_str() );
    if( !fin.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file " + filename ) );

    std::string str;

    int linec = 0;
    while( !fin.eof() ) {

	// Read line
	std::getline( fin, str );
	linec++;

	// Parse line
	if( str[0] == '#' )
	    continue;
	if( !parse_line( str, c, xscale, fscale, cdim, fdim, filename, linec ) )
	    continue;
	rec++;

	if( rec == 1 ) // Save first record
	    memcpy( c2, c, 6*sizeof(double) );
	else if( rec == 2 ) { // Calculate h 
	    for( i = 0; i < cdim; i++ ) {
		if( c[i] != c2[i] ) {
		    h = fabs(c2[i]-c[i]);
		    break;
		}
	    }
	    if( i == cdim )
		throw( Error( ERROR_LOCATION, "invalid coordinate data in " + filename + 
			      ":" + to_string(linec) + ": no difference in coordinate data" ) );
	}

	// Find origo (minimum coordinate values) and max (maximum coordinate values)
	for( i = 0; i < cdim; i++ ) {
	    if( c[i] < origo[i] )
		origo[i] = c[i];
	    if( c[i] > max[i] )
		max[i] = c[i];
	}
    }

    // Calculate size
    for( i = 0; i < cdim; i++ ) {
	size[i] = (size_t)floor((max[i]-origo[i])/h+0.5)+1;
	//if( size[i] == 1 )
	//throw( Error( ERROR_LOCATION, "vector field doesn\'t support maps with size 1" ) );
    }
    for( ; i < 3; i++ ) {
	size[i] = 1;
	origo[i] = 0.0;
	max[i] = 0.0;
    }

    // Check number of records
    if( rec != size[0]*size[1]*size[2] ) {
	throw( Error( ERROR_LOCATION, "number of records " + to_string(rec) + " in file " + filename + 
		      " doesn\'t match expected mesh size " + to_string(size[0]) + "x"
		      + to_string(size[1]) + "x" + to_string(size[2]) ) );
    }

    if( ibsimu.get_verbose_output() ) {
	Int3D one(1,1,1);
	std::cout << "  origo = " << origo << "\n";
	std::cout << "  size  = " << size << "\n";
	std::cout << "  max   = " << origo+h*(size-one) << "\n";
	std::cout << "  h     = " << h << "\n";
    }

    // Prepare VectorField
    _geom_mode = geom_mode;
    _size = size;
    _origo = origo;
    _h = h;
    _div_h = 1.0/_h;
    for( i = 0; i < 3; i++ ) {
	if( fout[i] )
	    _F[i] = new double[_size[0]*_size[1]*_size[2]];
	else
	    _F[i] = NULL;
     }

    // Read data to fill mesh
    int ind[3] = {0, 0, 0};
    fin.clear();
    fin.seekg( 0 );

    linec = 0;
    while( !fin.eof() ) {

	// Read line
	std::getline( fin, str );
	linec++;

	// Parse line
	if( str[0] == '#' )
	    continue;
	if( !parse_line( str, c, xscale, fscale, cdim, fdim, filename, linec ) )
	    continue;

	// Convert coordinates to indexes
	for( i = 0; i < cdim; i++ )
	    ind[i] = (int)floor((c[i]-origo[i])/h+0.5);

	size_t j;
	for( i = 0, j = 0; i < 3; i++ )
	    if( fout[i] ) {
		_F[i][size[0]*size[1]*ind[2] + size[0]*ind[1] + ind[0]] = c[3+j];
		j++;
	    }
    }

    fin.close();
}


VectorField::VectorField( const VectorField &f )
    : _geom_mode(f._geom_mode), _size(f._size), _origo(f._origo), _h(f._h), 
      _div_h(f._div_h)
{
    for( size_t i = 0; i < 3; i++ ) {
	if( f._F[i] != NULL ) {
	    _F[i] = new double[_size[0]*_size[1]*_size[2]];
	    memcpy( _F[i], f._F[i], _size[0]*_size[1]*_size[2]*sizeof(double) );
	} else {
	    _F[i] = NULL;
	}
    }
}


VectorField::~VectorField()
{
    for( size_t i = 0; i < 3; i++ ) {
	if( _F[i] != NULL )
	    delete [] _F[i];
    }
}


void VectorField::translate( Vec3D x )
{
    // Translate origo
    _origo += x;
}


void VectorField::transform( int ind[3] )
{
    Vec3D norigo;
    Int3D nsize;
    int p[3] = { abs(-ind[0])-1,
		 abs(-ind[1])-1,
		 abs(-ind[2])-1 };

    // Transform origo and size
    double *t[3];
    for( int32_t a = 0; a < 3; a++ ) {
	if( ind[a] < 0 )
	    norigo[a] = -_origo[p[a]] - (_size[p[a]]-1)*_h;
	else
	    norigo[a] = _origo[p[a]];
	nsize[a] = _size[p[a]];
    }

    // Transform mesh
    for( int32_t a = 0; a < 3; a++ ) {
	if( _F[p[a]] == NULL ) {
	    t[a] = NULL;
	    continue;
	} else {
	    // Allocate space for new mesh
	    t[a] = new double[_size[0]*_size[1]*_size[2]];
	}

	int32_t ni[3];
	int32_t i[3];
	double sign = ind[a] < 0 ? -1.0 : 1.0;
	for( ni[0] = 0; ni[0] < nsize[0]; ni[0]++ ) {
	    if( ind[0] < 0 )
		i[p[0]] = nsize[0]-1-ni[0];
	    else
		i[p[0]] = ni[0];

	    for( ni[1] = 0; ni[1] < nsize[1]; ni[1]++ ) {
		if( ind[1] < 0 )
		    i[p[1]] = nsize[1]-1-ni[1];
		else
		    i[p[1]] = ni[1];

		for( ni[2] = 0; ni[2] < nsize[2]; ni[2]++ ) {
		    if( ind[2] < 0 )
			i[p[2]] = nsize[2]-1-ni[2];
		    else
			i[p[2]] = ni[2];

		    t[a][(nsize[1]*ni[2] + ni[1])*nsize[0] + ni[0]] =
			sign*_F[p[a]][(_size[1]*i[2] + i[1])*_size[0] + i[0]];
		}
	    }
	}
    }

    // Replace field meshes with new ones
    for( int32_t a = 0; a < 3; a++ ) {
	if( _F[a] != NULL )
	    delete [] _F[a];
	_F[a] = t[a];
    }

    // Set origo and size
    _size = nsize;
    _origo = norigo;
}


void VectorField::scale( double s )
{
    if( s > 0 ) {
	_origo *= s;
	_h *= s;
	_div_h = 1.0/_h;
    } else if( s < 0 ) {
	int ind[3] = {-1, -2, -3};
	transform( ind );
	_origo *= fabs(s);
	_h *= fabs(s);
	_div_h = 1.0/_h;
    } else {
	throw( Error( ERROR_LOCATION, "scaling with zero" ) );
    }
}


void VectorField::rotate_x( int a )
{
    int r = a%360;
    if( r < 0 )
	r += 360;

    if( r == 0 ) {
	return;
    } else if( r == 90 ) {
	int ind[3] = {+1, -3, +2};
	transform( ind );
    } else if( r == 180 ) {
	int ind[3] = {+1, -2, -3};
	transform( ind );
    } else if( r == 270 ) {
	int ind[3] = {+1, +3, -2};
	transform( ind );
    } else {
	throw( Error( ERROR_LOCATION, "rotation angle not a multiple of 90" ) );
    }
}


void VectorField::rotate_y( int a )
{
    int r = a%360;
    if( r < 0 )
	r += 360;

    if( r == 0 ) {
	return;
    } else if( r == 90 ) {
	int ind[3] = {+3, +2, -1};
	transform( ind );
    } else if( r == 180 ) {
	int ind[3] = {-1, +2, -3};
	transform( ind );
    } else if( r == 270 ) {
	int ind[3] = {-3, +2, +1};
	transform( ind );
    } else {
	throw( Error( ERROR_LOCATION, "rotation angle not a multiple of 90.0" ) );
    }
}


void VectorField::rotate_z( int a )
{
    int r = a%360;
    if( r < 0 )
	r += 360;

    if( r == 0 ) {
	return;
    } else if( r == 90 ) {
	int ind[3] = {-2, +1, +3};
	transform( ind );
    } else if( r == 180 ) {
	int ind[3] = {-1, -2, +3};
	transform( ind );
    } else if( r == 270 ) {
	int ind[3] = {+2, -1, +3};
	transform( ind );
    } else {
	throw( Error( ERROR_LOCATION, "rotation angle not a multiple of 90.0" ) );
    }
}


void VectorField::clear()
{
    for( size_t i = 0; i < 3; i++ ) {
	if( _F[i] != NULL )
	    memset( _F[i], 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
    }
}


void VectorField::reset( geom_mode_e geom_mode, bool fout[3], Int3D size, 
			 Vec3D origo, double h )
{
    _geom_mode = geom_mode;
    _size      = size;
    _origo     = origo;
    _h         = h;
    _div_h     = 1.0/h;

    // Check mesh size legality
    if( _size[0] < 1 || _size[1] < 1 || _size[2] < 1 )
	throw( Error( ERROR_LOCATION, "illegal mesh size" ) );

    for( size_t i = 0; i < 3; i++ ) {
	if( _F[i] != NULL )
	    delete [] _F[i];
	if( fout[i] ) {
	    _F[i] = new double[_size[0]*_size[1]*_size[2]];
	    memset( _F[i], 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
	} else {
	    _F[i] = NULL;
	}
    }
}


void VectorField::get_minmax( double &min, double &max ) const
{
    min = std::numeric_limits<double>::infinity();
    max = -std::numeric_limits<double>::infinity();
    double val;

    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t a = 1; a < ncount; a++ ) {
	val = (*this)( a ).ssqr();
	if( val < min )
	    min = val;
	if( val > max )
	    max = val;
    }
    min = sqrt( min );
    max = sqrt( max );
}


void VectorField::get_defined_components( bool fout[3] ) const
{
    for( size_t a = 0; a < 3; a++ ) {
	if( _F[a] == NULL )
	    fout[a] = false;
	else
	    fout[a] = true;
    }
}


VectorField &VectorField::operator=( const VectorField &f )
{
    _geom_mode = f._geom_mode;
    _size      = f._size;
    _origo     = f._origo;
    _h         = f._h;
    _div_h     = f._div_h;
    for( size_t i = 0; i < 3; i++ ) {
	if( _F[i] != NULL )
	    delete [] _F[i];
	if( f._F[i] != NULL ) {
	    _F[i] = new double[_size[0]*_size[1]*_size[2]];
	    memcpy( _F[i], f._F[i], _size[0]*_size[1]*_size[2]*sizeof(double) );
	} else {
	    _F[i] = NULL;
	}
    }
    return( *this );
}


VectorField &VectorField::operator+=( const VectorField &f )
{
    if( _geom_mode != f._geom_mode || _size != f._size ||
	_origo != f._origo || _h != f._h )
	throw( Error( ERROR_LOCATION, "non-matching fields" ) );
    for( size_t i = 0; i < 3; i++ ) {
	if( (_F[i] == NULL) != (f._F[i] == NULL) )
	    throw( Error( ERROR_LOCATION, "non-matching fields" ) );
    }
    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t i = 0; i < 3; i++ ) {
	if( _F[i] != NULL ) {
	    for( size_t a = 0; a < ncount; a++ )
		_F[i][a] += f._F[i][a];
	}
    }
    return( *this );
}


VectorField &VectorField::operator*=( double x )
{
    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t i = 0; i < 3; i++ ) {
	if( _F[i] != NULL ) {
	    for( size_t a = 0; a < ncount; a++ )
		_F[i][a] *= x;
	}
    }
    return( *this );
}


VectorField &VectorField::operator/=( double x )
{
    double xi = 1.0/x;
    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t i = 0; i < 3; i++ ) {
	if( _F[i] != NULL ) {
	    for( size_t a = 0; a < ncount; a++ )
		_F[i][a] *= xi;
	}
    }
    return( *this );
}


const Vec3D VectorField::operator()( int32_t i ) const
{
    Vec3D ret;

    for( size_t b = 0; b < 3; b++ ) {
	if( _F[b] != NULL ) {
	    ret[b] = _F[b][i];
	}
    }

    return( ret );
}


const Vec3D VectorField::operator()( int32_t i, int32_t j ) const
{
    Vec3D ret;

    for( size_t b = 0; b < 3; b++ ) {
	if( _F[b] != NULL ) {
	    ret[b] = _F[b][i+j*_size[0]];
	}
    }

    return( ret );
}


const Vec3D VectorField::operator()( int32_t i, int32_t j, int32_t k ) const
{
    Vec3D ret;

    for( size_t b = 0; b < 3; b++ ) {
	if( _F[b] != NULL ) {
	    ret[b] = _F[b][i+(j+k*_size[1])*_size[0]];
	}
    }

    return( ret );
}


void VectorField::set( int32_t i, const Vec3D &v )
{
    for( size_t b = 0; b < 3; b++ ) {
	if( _F[b] != NULL ) {
	    _F[b][i] = v[b];
	}
    }
}
    

void VectorField::set( int32_t i, int32_t j, const Vec3D &v )
{
    for( size_t b = 0; b < 3; b++ ) {
	if( _F[b] != NULL ) {
	    _F[b][i+j*_size[0]] = v[b];
	}
    }
}


void VectorField::set( int32_t i, int32_t j, int32_t k, const Vec3D &v )
{
    for( size_t b = 0; b < 3; b++ ) {
	if( _F[b] != NULL ) {
	    _F[b][i+(j+k*_size[1])*_size[0]] = v[b];
	}
    }
}


Vec3D VectorField::operator()( Vec3D x ) const
{
    Vec3D ret;

    switch( _geom_mode ) {
    case MODE_1D:
    {
	// Constant field
	if( _size[0] == 1 ) {
	    for( size_t b = 0; b < 3; b++ ) {
		if( _F[b] != NULL ) {
		    ret[b] = _F[b][0];
		}
	    }
	    break;
	}

	// Limit to double the simulation box -> return zero
	if( x[0] < _origo[0]-_size[0]*_h ) {
	    break;
	} else if( x[0] > _origo[0]+2.0*_size[0]*_h ) {
	    break;
	}

	// Linear approximation
	int32_t i = (int32_t)floor( (x[0]-_origo[0])*_div_h );
	if( i < 0 )
	    i = 0;
	else if( i >= _size[0]-1 )
	    i = _size[0]-2;

	double t = _div_h*( x[0]-(i*_h+_origo[0]) );
	
	for( size_t b = 0; b < 3; b++ ) {
	    if( _F[b] != NULL ) {
		ret[b] = (1.0-t)*_F[b][i] + t*_F[b][i+1];
	    }
	}
	break;
    }
    case MODE_2D:
    case MODE_CYL:
    {
	int32_t i, j, di, dj;
	double t, u;

	// Limit to double the simulation box -> return zero
	if( x[0] < _origo[0]-_size[0]*_h ) {
	    break;
	} else if( x[0] > _origo[0]+2.0*_size[0]*_h ) {
	    break;
	} else if( x[1] < _origo[1]-_size[1]*_h ) {
	    break;
	} else if( x[1] > _origo[1]+2.0*_size[1]*_h ) {
	    break;
	}

	if( _size[0] == 1 ) {
	    i  = 0;
	    di = 0;
	    t  = 0.0;
	} else {
	    i = (int32_t)floor( (x[0]-_origo[0])*_div_h );
	    if( i < 0 )
		i = 0;
	    else if( i >= _size[0]-1 )
		i = _size[0]-2;
	    t  = _div_h*( x[0]-(i*_h+_origo[0]) );
	    di = 1;
	}

	if( _size[1] == 1 ) {
	    j  = 0;
	    dj = 0;
	    u  = 0.0;
	} else {
	    j = (int32_t)floor( (x[1]-_origo[1])*_div_h );
	    if( j < 0 )
		j = 0;
	    else if( j >= _size[1]-1 )
		j = _size[1]-2;
	    u = _div_h*( x[1]-(j*_h+_origo[1]) );
	    dj = _size[0];
	}

	int32_t base = _size[0]*j + i;
	for( size_t b = 0; b < 3; b++ ) {
	    if( _F[b] != NULL ) {
		ret[b] = (1.0-t)*(1.0-u)*_F[b][base] +
		         (    t)*(1.0-u)*_F[b][base+di] +
		         (1.0-t)*(    u)*_F[b][base+dj] +
		         (    t)*(    u)*_F[b][base+di+dj];
	    }
	}
	break;
    }
    default:
    {
	int32_t i, j, k, di, dj, dk;
	double t, u, v;

	// Limit to double the simulation box -> return zero
	if( x[0] < _origo[0]-_size[0]*_h ) {
	    break;
	} else if( x[0] > _origo[0]+2.0*_size[0]*_h ) {
	    break;
	} else if( x[1] < _origo[1]-_size[1]*_h ) {
	    break;
	} else if( x[1] > _origo[1]+2.0*_size[1]*_h ) {
	    break;
	} else if( x[2] < _origo[2]-_size[2]*_h ) {
	    break;
	} else if( x[2] > _origo[2]+2.0*_size[2]*_h ) {
	    break;
	}

	if( _size[0] == 1 ) {
	    i  = 0;
	    di = 0;
	    t  = 0.0;
	} else {
	    i = (int32_t)floor( (x[0]-_origo[0])*_div_h );
	    if( i < 0 )
		i = 0;
	    else if( i >= _size[0]-1 )
		i = _size[0]-2;
	    t  = _div_h*( x[0]-(i*_h+_origo[0]) );
	    di = 1;
	}

	if( _size[1] == 1 ) {
	    j  = 0;
	    dj = 0;
	    u  = 0.0;
	} else {
	    j = (int32_t)floor( (x[1]-_origo[1])*_div_h );
	    if( j < 0 )
		j = 0;
	    else if( j >= _size[1]-1 )
		j = _size[1]-2;
	    u = _div_h*( x[1]-(j*_h+_origo[1]) );
	    dj = _size[0];
	}

	if( _size[2] == 1 ) {
	    k  = 0;
	    dk = 0;
	    v  = 0.0;
	} else {
	    k = (int32_t)floor( (x[2]-_origo[2])*_div_h );
	    if( k < 0 )
		k = 0;
	    else if( k >= _size[2]-1 )
		k = _size[2]-2;
	    v = _div_h*( x[2]-(k*_h+_origo[2]) );
	    dk = _size[0]*_size[1];
	}

	int32_t base = (k*_size[1] + j)*_size[0] + i;
	for( size_t b = 0; b < 3; b++ ) {
	    if( _F[b] != NULL ) {
		ret[b] = (1.0-t)*(1.0-u)*(1.0-v)*_F[b][base] +
		         (    t)*(1.0-u)*(1.0-v)*_F[b][base+di] +
		         (1.0-t)*(    u)*(1.0-v)*_F[b][base+dj] +
		         (    t)*(    u)*(1.0-v)*_F[b][base+di+dj] +
		         (1.0-t)*(1.0-u)*(    v)*_F[b][base+dk] +
		         (    t)*(1.0-u)*(    v)*_F[b][base+di+dk] +
		         (1.0-t)*(    u)*(    v)*_F[b][base+dj+dk] +
		         (    t)*(    u)*(    v)*_F[b][base+di+dj+dk];
	    }
	}
	break;
    }
    }

    return( ret );
}


void VectorField::debug_print( void ) const
{
    int32_t a;
    std::cout << "**VectorField\n";
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
    for( size_t i = 0; i < 3; i++ ) {
	std::cout << "F[" << i << "] = ";
	if( _F[i] == NULL ) {
	    std::cout << "0\n";
	    continue;
	}
	std::cout << "(";
	for( a = 0; a < _size[0]*_size[1]*_size[2]-1; a++ )
	    std::cout << _F[i][a] << ", ";
	if( a < _size[0]*_size[1]*_size[2] )
	    std::cout << _F[i][a] << ")\n";
    }

}


















