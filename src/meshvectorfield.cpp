/*! \file meshvectorfield.cpp
 *  \brief %Mesh based vector fields
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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
#include <limits>
#include <cstdlib>
#include "readascii.hpp"
#include "meshvectorfield.hpp"
#include "ibsimu.hpp"


MeshVectorField::MeshVectorField()
{
    _F[0] = _F[1] = _F[2] = NULL;
    _extrpl[0] = _extrpl[1] = _extrpl[2] = _extrpl[3] = _extrpl[4] = _extrpl[5] = FIELD_EXTRAPOLATE;
}


MeshVectorField::MeshVectorField( const Mesh &m, const bool fout[3] )
    : Mesh(m)
{
    _extrpl[0] = _extrpl[1] = _extrpl[2] = _extrpl[3] = _extrpl[4] = _extrpl[5] = FIELD_EXTRAPOLATE;

    check_definition();

    for( size_t i = 0; i < 3; i++ ) {
	if( fout[i] ) {
	    _F[i] = new double[_size[0]*_size[1]*_size[2]];
	    memset( _F[i], 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
	} else {
	    _F[i] = NULL;
	}
    }
}


MeshVectorField::MeshVectorField( geom_mode_e geom_mode, const bool fout[3], Int3D size, 
				  Vec3D origo, double h )
    : Mesh(geom_mode,size,origo,h)
{
    _extrpl[0] = _extrpl[1] = _extrpl[2] = _extrpl[3] = _extrpl[4] = _extrpl[5] = FIELD_EXTRAPOLATE;

    check_definition();

    for( size_t i = 0; i < 3; i++ ) {
	if( fout[i] ) {
	    _F[i] = new double[_size[0]*_size[1]*_size[2]];
	    memset( _F[i], 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
	} else {
	    _F[i] = NULL;
	}
    }
}


MeshVectorField::MeshVectorField( std::istream &s )
    : Mesh(s)
{
    check_definition();

    for( int i = 0; i < 6; i++ )
	_extrpl[i] = (field_extrpl_e)read_int32( s );

    for( int i = 0; i < 3; i++ ) {
	if( read_int8( s ) ) {
	    _F[i] = new double[_size[0]*_size[1]*_size[2]];
	    read_compressed_block( s, _size[0]*_size[1]*_size[2]*sizeof(double), (int8_t *)_F[i] );
	} else {
	    _F[i] = NULL;
	}
    }
}


void MeshVectorField::check_definition()
{
    // Check mesh size legality
    if( _size[0] < 1 || _size[1] < 1 || _size[2] < 1 )
	throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
}


MeshVectorField::MeshVectorField( geom_mode_e geom_mode, const bool fout[3], double xscale, 
				  double fscale, const std::string &filename )
{
    _extrpl[0] = _extrpl[1] = _extrpl[2] = _extrpl[3] = _extrpl[4] = _extrpl[5] = FIELD_EXTRAPOLATE;

    if( ibsimu.get_verbose_output() )
	std::cout << "Reading vector field from " << filename << "\n";

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
    }

    // Set number of field components (fdim)
    int fdim = 0;
    for( int a = 0; a < 3; a++ ) {
	if( fout[a] )
	    fdim++;
    }
    if( fdim == 0 )
	throw( Error( ERROR_LOCATION, "no field components to read" ) );

    // Read data columns
    int columns = cdim + fdim;
    ReadAscii data( filename, columns );

    // Check size
    if( data.rows() < 2 )
	throw( Error( ERROR_LOCATION, "too few records read" ) );
    
    // Determine mesh parameters from data
    Int3D size;
    Vec3D origo( std::numeric_limits<double>::infinity(), 
		 std::numeric_limits<double>::infinity(), 
		 std::numeric_limits<double>::infinity() );
    Vec3D max( -std::numeric_limits<double>::infinity(), 
	       -std::numeric_limits<double>::infinity(), 
	       -std::numeric_limits<double>::infinity() );
    
    // Find origo (minimum coordinate values) and max (maximum coordinate values)
    for( uint32_t a = 0; a < data.rows(); a++ ) {
	for( uint32_t b = 0; b < cdim; b++ ) {
	    double d = xscale*data[b][a];
	    if( d < origo[b] )
		origo[b] = d;
	    if( d > max[b] )
		max[b] = d;
	}
    }
    
    // Initial guess for h is difference of first two records
    uint32_t first_dir = 0;
    double h = 0.0;
    for( uint32_t b = 0; b < cdim; b++ ) {
	if( fabs(data[b][1] - data[b][0]) != 0.0 ) {
	    first_dir = b;
	    h = xscale*fabs(data[b][1] - data[b][0]);
	    break;
	}
    }
    if( first_dir == cdim )
	throw( Error( ERROR_LOCATION, "no coordinate step recognized" ) );	
    
    // Calculate size from h
    uint32_t i;
    for( i = 0; i < cdim; i++ )
	size[i] = (size_t)floor((max[i]-origo[i])/h+0.5)+1;
    for( ; i < 3; i++ ) {
	size[i] = 1;
	origo[i] = 0.0;
	max[i] = 0.0;
    }

    // Correct h to higher accuracy
    h = (max[first_dir]-origo[first_dir]) / ((double)size[first_dir]-1.0);

    // Check number of records
    if( data.rows() != (uint32_t)size[0]*size[1]*size[2] ) {
	throw( Error( ERROR_LOCATION, "number of records " + to_string(data.rows()) +
		      " in file " + filename + " doesn\'t match expected mesh size " +
		      to_string(size[0]) + "x" + to_string(size[1]) + "x" + to_string(size[2]) ) );
    }

    if( ibsimu.get_verbose_output() ) {
	std::cout << "  origo = " << origo << "\n";
	std::cout << "  size  = " << size << "\n";
	std::cout << "  max   = " << max << "\n";
	std::cout << "  h     = " << h << "\n";
    }

    // Prepare MeshVectorField
    Mesh::reset( geom_mode, size, origo, h );
    check_definition();
    for( i = 0; i < 3; i++ ) {
	if( fout[i] )
	    _F[i] = new double[_size[0]*_size[1]*_size[2]];
	else
	    _F[i] = NULL;
    }

    // Read data to fill mesh
    int ind[3] = {0, 0, 0};
    for( uint32_t a = 0; a < data.rows(); a++ ) {

	// Convert coordinates to indexes
	for( uint32_t b = 0; b < cdim; b++ )
	    ind[b] = (int)floor((xscale*data[b][a]-origo[b])/h+0.5);

	uint32_t j = 0;
	for( uint32_t i = 0; i < 3; i++ ) {
	    if( fout[i] ) {
		_F[i][size[0]*size[1]*ind[2] + size[0]*ind[1] + ind[0]] = fscale*data[cdim+j][a];
		j++;
	    }
	}
    }
}


MeshVectorField::MeshVectorField( const MeshVectorField &f )
    : Mesh(f)
{
    for( size_t i = 0; i < 6; i++ )
	_extrpl[i] = f._extrpl[i];

    for( size_t i = 0; i < 3; i++ ) {
	if( f._F[i] != NULL ) {
	    _F[i] = new double[_size[0]*_size[1]*_size[2]];
	    memcpy( _F[i], f._F[i], _size[0]*_size[1]*_size[2]*sizeof(double) );
	} else {
	    _F[i] = NULL;
	}
    }
}


MeshVectorField::~MeshVectorField()
{
    for( size_t i = 0; i < 3; i++ ) {
	if( _F[i] != NULL )
	    delete [] _F[i];
    }
}


void MeshVectorField::translate( Vec3D x )
{
    // Translate origo
    _origo += x;

    // Translate max
    _max += x;
}


void MeshVectorField::transform( int ind[3] )
{
    Vec3D norigo;
    Vec3D nmax;
    Int3D nsize;
    int p[3] = { abs(-ind[0])-1,
		 abs(-ind[1])-1,
		 abs(-ind[2])-1 };

    // Transform origo, size and max
    double *t[3];
    for( int32_t a = 0; a < 3; a++ ) {
	if( ind[a] < 0 )
	    norigo[a] = -_origo[p[a]] - (_size[p[a]]-1)*_h;
	else
	    norigo[a] = _origo[p[a]];
	nsize[a] = _size[p[a]];
	nmax[a] = _origo[a]+_h*(_size[a]-1);
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
    _max = nmax;
    _origo = norigo;
}


void MeshVectorField::scale( double s )
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


void MeshVectorField::rotate_x( int a )
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


void MeshVectorField::rotate_y( int a )
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


void MeshVectorField::rotate_z( int a )
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


void MeshVectorField::clear()
{
    for( size_t i = 0; i < 3; i++ ) {
	if( _F[i] != NULL )
	    memset( _F[i], 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
    }
}


void MeshVectorField::reset( geom_mode_e geom_mode, const bool fout[3], Int3D size, 
			 Vec3D origo, double h )
{
    Mesh::reset( geom_mode, size, origo, h );
    check_definition();

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


void MeshVectorField::get_minmax( double &min, double &max ) const
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


void MeshVectorField::get_defined_components( bool fout[3] ) const
{
    for( size_t a = 0; a < 3; a++ ) {
	if( _F[a] == NULL )
	    fout[a] = false;
	else
	    fout[a] = true;
    }
}


MeshVectorField &MeshVectorField::operator=( const MeshVectorField &f )
{
    (Mesh)(*this) = (Mesh)f;

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


MeshVectorField &MeshVectorField::operator+=( const MeshVectorField &f )
{
    if( (Mesh)(*this) == (Mesh)f )
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


MeshVectorField &MeshVectorField::operator*=( double x )
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


MeshVectorField &MeshVectorField::operator/=( double x )
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


const Vec3D MeshVectorField::operator()( int32_t i ) const
{
    Vec3D ret;

    for( size_t b = 0; b < 3; b++ ) {
	if( _F[b] != NULL ) {
	    ret[b] = _F[b][i];
	}
    }

    return( ret );
}


const Vec3D MeshVectorField::operator()( int32_t i, int32_t j ) const
{
    Vec3D ret;

    for( size_t b = 0; b < 3; b++ ) {
	if( _F[b] != NULL ) {
	    ret[b] = _F[b][i+j*_size[0]];
	}
    }

    return( ret );
}


const Vec3D MeshVectorField::operator()( int32_t i, int32_t j, int32_t k ) const
{
    Vec3D ret;

    for( size_t b = 0; b < 3; b++ ) {
	if( _F[b] != NULL ) {
	    ret[b] = _F[b][i+(j+k*_size[1])*_size[0]];
	}
    }

    return( ret );
}


void MeshVectorField::set( int32_t i, const Vec3D &v )
{
    for( size_t b = 0; b < 3; b++ ) {
	if( _F[b] != NULL ) {
	    _F[b][i] = v[b];
	}
    }
}
    

void MeshVectorField::set( int32_t i, int32_t j, const Vec3D &v )
{
    for( size_t b = 0; b < 3; b++ ) {
	if( _F[b] != NULL ) {
	    _F[b][i+j*_size[0]] = v[b];
	}
    }
}


void MeshVectorField::set( int32_t i, int32_t j, int32_t k, const Vec3D &v )
{
    for( size_t b = 0; b < 3; b++ ) {
	if( _F[b] != NULL ) {
	    _F[b][i+(j+k*_size[1])*_size[0]] = v[b];
	}
    }
}


const Vec3D MeshVectorField::operator()( Vec3D x ) const
{
    Vec3D R;
    Vec3D sign( 1.0, 1.0, 1.0 );

    if( _size[0] == 0 )
	return( R );

    switch( _geom_mode ) {
    case MODE_1D:
    {
	// Constant field
	if( _size[0] == 1 ) {
	    for( size_t b = 0; b < 3; b++ ) {
		if( _F[b] != NULL ) {
		    R[b] = _F[b][0];
		}
	    }
	    break;
	}

	if( x[0] < _origo[0] ) {
	    if( _extrpl[0] == FIELD_ZERO ) {
		// return zero
		return( Vec3D(0.0) );
	    } else if( _extrpl[0] == FIELD_NAN ) {
		// return NaN
		return( Vec3D(std::numeric_limits<double>::quiet_NaN()) );
	    }
	    if( x[0] < _origo[0]-_size[0]*_h ) {
		// Outside double the simulation box: return zero
		return( Vec3D(0.0) );
	    }
	    if( _extrpl[0] == FIELD_MIRROR ) {
		sign[0] *= -1.0;
		x[0]  = 2.0*_origo[0] - x[0];
	    }
		
	} else if( x[0] > _max[0] ) {
	    if( _extrpl[1] == FIELD_ZERO ) {
		// return zero
		return( Vec3D(0.0) );
	    } else if( _extrpl[1] == FIELD_NAN ) {
		// return NaN
		return( Vec3D(std::numeric_limits<double>::quiet_NaN()) );
	    }
	    if( x[0] > _origo[0]+2.0*_size[0]*_h ) {
		// Outside double the simulation box: return zero
		return( R );
	    }
	    if( _extrpl[1] == FIELD_MIRROR ) {
		sign[0] *= -1.0;
		x[0]  = 2.0*_max[0] - x[0];
	    }
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
		R[b] = sign[b]*( (1.0-t)*_F[b][i] + t*_F[b][i+1] );
	    }
	}
	break;
    }
    case MODE_2D:
    case MODE_CYL:
    {
	int32_t i, j, di, dj;
	double t, u;

	for( int a = 0; a < 2; a++ ) {
	    if( x[a] < _origo[a] ) {
		if( _extrpl[2*a] == FIELD_ZERO ) {
                    // return zero
                    return( Vec3D(0.0) );
                } else if( _extrpl[2*a] == FIELD_NAN ) {
		    // return NaN
		    return( Vec3D(std::numeric_limits<double>::quiet_NaN()) );
		}
		if( x[a] < _origo[a]-_size[a]*_h ) {
		    // Limit to double the simulation box -> return zero
		    return( R );
		}
		if( _extrpl[2*a] == FIELD_MIRROR ) {
                    sign[a] *= -1.0;
                    x[a]     = 2.0*_origo[a] - x[a];
		}

	    } else if( x[0] > _max[0] ) {
		if( _extrpl[2*a+1] == FIELD_ZERO ) {
                    // return zero
                    return( Vec3D(0.0) );
                } else if( _extrpl[2*a+1] == FIELD_NAN ) {
		    // return NaN
		    return( Vec3D(std::numeric_limits<double>::quiet_NaN()) );
		}
		if( x[a] > _origo[a]+2.0*_size[a]*_h ) {
		    // Limit to double the simulation box -> return zero
		    return( R );
		}
                if( _extrpl[2*a+1] == FIELD_MIRROR ) {
                    sign[a] *= -1.0;
                    x[a]     = 2.0*_max[a] - x[a];
		}
	    }
	}

	if( _size[0] == 1 ) {
	    // Field constant in x-direction
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
	    // Field constant in y-direction
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
		R[b] = sign[b]*( (1.0-t)*(1.0-u)*_F[b][base] +
				 (    t)*(1.0-u)*_F[b][base+di] +
				 (1.0-t)*(    u)*_F[b][base+dj] +
				 (    t)*(    u)*_F[b][base+di+dj] );
	    }
	}
	break;
    }
    default:
    {
	int32_t i, j, k, di, dj, dk;
	double t, u, v;

	for( int a = 0; a < 3; a++ ) {
	    if( x[a] < _origo[a] ) {
		if( _extrpl[2*a] == FIELD_ZERO ) {
                    // return zero
                    return( Vec3D(0.0) );
                } else if( _extrpl[2*a] == FIELD_NAN ) {
		    // return NaN
		    return( Vec3D(std::numeric_limits<double>::quiet_NaN()) );
		}
		if( x[a] < _origo[a]-_size[a]*_h ) {
		    // Limit to double the simulation box -> return zero
		    return( R );
		}
		if( _extrpl[2*a] == FIELD_MIRROR ) {
                    sign[a] *= -1.0;
                    x[a]     = 2.0*_origo[a] - x[a];
		}

	    } else if( x[0] > _max[0] ) {
		if( _extrpl[2*a+1] == FIELD_ZERO ) {
                    // return zero
                    return( Vec3D(0.0) );
                } else if( _extrpl[2*a+1] == FIELD_NAN ) {
		    // return NaN
		    return( Vec3D(std::numeric_limits<double>::quiet_NaN()) );
		}
		if( x[a] > _origo[a]+2.0*_size[a]*_h ) {
		    // Limit to double the simulation box -> return zero
		    return( R );
		}
                if( _extrpl[2*a+1] == FIELD_MIRROR ) {
                    sign[a] *= -1.0;
                    x[a]     = 2.0*_max[a] - x[a];
		}
            }
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
		R[b] = sign[b]*( (1.0-t)*(1.0-u)*(1.0-v)*_F[b][base] +
				 (    t)*(1.0-u)*(1.0-v)*_F[b][base+di] +
				 (1.0-t)*(    u)*(1.0-v)*_F[b][base+dj] +
				 (    t)*(    u)*(1.0-v)*_F[b][base+di+dj] +
				 (1.0-t)*(1.0-u)*(    v)*_F[b][base+dk] +
 				 (    t)*(1.0-u)*(    v)*_F[b][base+di+dk] +
				 (1.0-t)*(    u)*(    v)*_F[b][base+dj+dk] +
				 (    t)*(    u)*(    v)*_F[b][base+di+dj+dk] );
	    }
	}
	break;
    }
    }

    return( R );
}


void MeshVectorField::save( std::ostream &s ) const
{
    Mesh::save( s );

    for( int i = 0; i < 6; i++ )
	write_int32( s, _extrpl[i] );

    for( int i = 0; i < 3; i++ ) {
	if( _F[i] ) {
	    write_int8( s, 1 );
	    write_compressed_block( s, _size[0]*_size[1]*_size[2]*sizeof(double), (int8_t *)_F );
	} else {
	    write_int8( s, 0 );
	}
    }
}


void MeshVectorField::debug_print( std::ostream &os ) const
{
    Mesh::debug_print( os );

    std::cout << "**MeshVectorField\n";
    std::cout << "extrpl = (" 
	      << _extrpl[0] << ", "
	      << _extrpl[1] << ", "
	      << _extrpl[2] << ", "
	      << _extrpl[3] << ", "
	      << _extrpl[4] << ", "
	      << _extrpl[5] << ")\n";
    for( size_t i = 0; i < 3; i++ ) {
	std::cout << "F[" << i << "] = ";
	if( _F[i] == NULL ) {
	    std::cout << "NULL\n";
	    continue;
	}
	std::cout << "(";
	if( _size[0]*_size[1]*_size[2] < 10 ) {
	    int a;
	    for( a = 0; a < _size[0]*_size[1]*_size[2]-1; a++ )
		std::cout << _F[i][a] << ", ";
	    if( a < _size[0]*_size[1]*_size[2] )
		std::cout << _F[i][a] << ")\n";
	} else {
	    // Print only 10 first nodes
	    for( int a = 0; a < 10; a++ )
		std::cout << _F[i][a] << ", ";
	    std::cout << "... )\n";
	}
    }
}

