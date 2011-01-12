/*! \file colormap.cpp
 *  \brief Source code for colormap.cpp
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

#include <cmath>
#include <limits>
#include "compmath.hpp"
#include "colormap.hpp"


Colormap::Colormap()
    : _interpolation(INTERPOLATION_BILINEAR), _zscale(ZSCALE_LINEAR), _n(0), _m(0), _intrp(NULL)
{

}


Colormap::Colormap( const double datarange[4], size_t n, size_t m, 
		    const std::vector<double> &data )
    : _interpolation(INTERPOLATION_BILINEAR), _zscale(ZSCALE_LINEAR), _n(n), _m(m), _intrp(NULL)
{
    _datarange[0] = datarange[0];
    _datarange[1] = datarange[1];
    _datarange[2] = datarange[2];
    _datarange[3] = datarange[3];
    
    if( n*m != data.size() )
	throw( Error( ERROR_LOCATION, "data size not equal to n*m" ) );
    _f = data;

    // Go through data to find minimum and maximum
    _zmin = std::numeric_limits<double>::infinity();
    _zmax = -std::numeric_limits<double>::infinity();
    for( size_t j = 0; j < _m; j++ ) {
	for( size_t i = 0; i < _n; i++ ) {
	    int p = i+j*_n;
	    if( _f[p] < _zmin )
		_zmin = _f[p];
	    if( _f[p] > _zmax )
		_zmax = _f[p];
	}
    }
    if( _zmin == _zmax ) {
	_zmax += 1.0;
	//_zmin -= 1.0;
    }

    make_data_interpolation();
}


Colormap::~Colormap()
{
    if( _intrp )
	delete _intrp;
}


void Colormap::set_interpolation( interpolation_e interpolation )
{
    _interpolation = interpolation;
    make_data_interpolation();
}


void Colormap::set_zscale( zscale_e zscale )
{
    _zscale = zscale;
}


void Colormap::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
    // Prepare by fetching surface and its parameters
    cairo_surface_t *surface = cairo_get_target( cairo );
    unsigned char *buf = cairo_image_surface_get_data( surface );
    cairo_format_t format = cairo_image_surface_get_format( surface );
    if( buf == NULL ) {
	throw( Error( ERROR_LOCATION, "cairo surface not image surface" ) );
    } else if( format != CAIRO_FORMAT_ARGB32 ) {
	throw( Error( ERROR_LOCATION, "cairo image surface format not supported" ) );
    }
    int width  = cairo_image_surface_get_width( surface );
    int height = cairo_image_surface_get_height( surface );
    int stride = cairo_image_surface_get_stride( surface );

    // Calculate pixel ranges
    double prange[4];
    if( _datarange[0] < range[0] )
	prange[0] = range[0];
    else
	prange[0] = _datarange[0];
    if( _datarange[1] < range[1] )
	prange[1] = range[1];
    else
	prange[1] = _datarange[1];

    if( _datarange[2] > range[2] )
	prange[2] = range[2];
    else
	prange[2] = _datarange[2];
    if( _datarange[3] > range[3] )
	prange[3] = range[3];
    else
	prange[3] = _datarange[3];

    /*
    std::cout << "Untransformed prange: "
	      << prange[0] << " "
	      << prange[1] << " "
	      << prange[2] << " "
	      << prange[3] << "\n";
    */

    cm->transform( prange[0], prange[1] );
    cm->transform( prange[2], prange[3] );

    /*
    std::cout << "Transformed prange: "
	      << prange[0] << " "
	      << prange[1] << " "
	      << prange[2] << " "
	      << prange[3] << "\n";
    */

    // Calculate pixel limits of drawn area, y flipped to have smaller
    // numbers as 0 and 1, bigger as 2 and 3.
    int plim[4] = { (int)floor(prange[0]+0.5),
		    (int)floor(prange[3]+0.5),
		    (int)floor(prange[2]+0.5),
		    (int)floor(prange[1]+0.5) };

    // Check that drawing is only done to valid buffer
    if( plim[0] < 0 )
	plim[0] = 0;
    if( plim[1] < 0 )
	plim[1] = 0;
    if( plim[2] >= width )
	plim[2] = width-1;
    if( plim[3] >= height )
	plim[3] = height-1;
    if( plim[0] > plim[2] || plim[1] > plim[3] )
	return;

    /*
    std::cout << width << "\n";
    std::cout << height << "\n";
    std::cout << stride << "\n";
    std::cout << "Rounded prange: " 
	      << plim[0] << " "
	      << plim[1] << " "
	      << plim[2] << " "
	      << plim[3] << "\n";
    */

    if( !_intrp )
	throw( Error( ERROR_LOCATION, "No data available" ) );

    // Go through pixel limits
    // Error if either end is at zero with LOG scaling
    if( _zscale == ZSCALE_LOG && _zmin <= 0.0 && _zmax >= 0.0 )
	throw( Error( ERROR_LOCATION, "zmin and zmax on different sides of zero" ) );

    int sign;
    if( _zmax <= 0.0 ) {
	// Completely on negative side
	sign = -1;
    } else if( _zmin >= 0.0 ) {
	// Completely on positive side
	sign = +1;
    } else {
	// Both negative and positive
	sign = 0;
    }
	
    for( int i = plim[0]; i <= plim[2]; i++ ) {
	for( int j = plim[1]; j <= plim[3]; j++ ) {

	    // Tranform to logical coordinates
	    double x[2] = { i, j };
	    cm->inv_transform( x[0], x[1] );

	    double val = get_value( x[0], x[1] );

	    // Scale value to [0:1] maintaining monotonic rising property
	    if( _zscale == ZSCALE_LINEAR )
		val = (val-_zmin)/(_zmax-_zmin);
	    else if( _zscale == ZSCALE_LOG ) {
		if( sign > 0 )
		    val = (log(val)-log(_zmin)) / (log(_zmax)-log(_zmin));
		else
		    val = (log(-val)-log(-_zmin)) / (log(-_zmax)-log(-_zmin));
	    } else if( _zscale == ZSCALE_RELLOG ) {
		if( sign > 0 ) {
		    val = (val-_zmin)/(_zmax-_zmin);
		    val = (log(0.001+val) - log(0.001)) / (log(1.001) - log(0.001));
		} else if( sign < 0 ) {
		    val = (val-_zmax)/(_zmin-_zmax);
		    val = 1.0+(log(0.001) - log(0.001+val)) / (log(1.001) - log(0.001));
		} else {
		    if( val > 0.0 ) {
			val = (val-0.0)/(_zmax-0.0);
			val = 0.5+0.5*(log(0.001+val) - log(0.001)) / (log(1.001) - log(0.001));
		    } else {
			val = (val-0.0)/(_zmin-0.0);
			val = 0.5+0.5*(log(0.001) - log(0.001+val)) / (log(1.001) - log(0.001));
		    }
		}
	    }
	    Color c;
	    if( comp_isinf( val ) || comp_isnan( val ) )
		continue;

	    c = _palette( val );
	    buf[j*stride+4*i+0] = (unsigned char)(255*c[2]);  // Blue
	    buf[j*stride+4*i+1] = (unsigned char)(255*c[1]);  // Green
	    buf[j*stride+4*i+2] = (unsigned char)(255*c[0]);  // Red
	    buf[j*stride+4*i+3] = (unsigned char)255;         // Alpha
	}
    }
}


void Colormap::get_bbox( double bbox[4] )
{
    bbox[0] = _datarange[0];
    bbox[1] = _datarange[1];
    bbox[2] = _datarange[2];
    bbox[3] = _datarange[3];
}


void Colormap::set_palette( const Palette &palette )
{
    _palette = palette;
}


void Colormap::get_zrange( double &min, double &max ) const
{
    min = _zmin;
    max = _zmax;
}


void Colormap::make_data_interpolation( void )
{
    // Free old interpolation
    if( _intrp )
	delete _intrp;

    // Make a new interpolation
    switch( _interpolation ) {
    case INTERPOLATION_CLOSEST:
	_intrp = new ClosestInterpolation2D( _n, _m, _f );
	break;
    case INTERPOLATION_BILINEAR:
	_intrp = new BiLinearInterpolation2D( _n, _m, _f );
	break;
    case INTERPOLATION_BICUBIC:
	_intrp = new BiCubicInterpolation2D( _n, _m, _f );
	break;
    default:
	throw( Error( ERROR_LOCATION, "unknown interpoltaion type" ) );	
    }
}


double Colormap::get_value( double x, double y ) const
{	    
    // Calculate relative point in data
    double t = (x-_datarange[0])/(_datarange[2]-_datarange[0]);
    double u = (y-_datarange[1])/(_datarange[3]-_datarange[1]);
    
    return( (*_intrp)( t, u ) );
}



