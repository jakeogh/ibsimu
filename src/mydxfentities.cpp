/*! \file mydxfentities.cpp
 *  \brief DXF entities
 */

/* Copyright (c) 2010 Taneli Kalvas. All rights reserved.
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


#include <iostream>
#include <math.h>
#include <limits>
#include "mydxfentities.hpp"
#include "polysolver.hpp"
#include "error.hpp"


#define MAX_ARC_CHOP 4

//#define MYDXF_DEBUG_BBOX



//#define MYDXF_DEBUG 1

/* ************************************************************************** *
 * DXFEntity                                                                  *
 * ************************************************************************** */


MyDXFEntity::MyDXFEntity()
{

}


void MyDXFEntity::bbox_ppoint( Vec3D &min, Vec3D &max, const Vec3D &p )
{
    for( int a = 0; a < 3; a++ ) {
	if( p[a] < min[a] )
	    min[a] = p[a];
	if( p[a] > max[a] )
	    max[a] = p[a];
    }
}


void MyDXFEntity::debug_print_base( std::ostream &os ) const
{
    os << "  handle = \'" << _handle << "\'\n";
    os << "  layer = \'" << _layer << "\'\n";
}


std::ostream &operator<<( std::ostream &os, const MyDXFEntity &ent )
{
    ent.debug_print( os );
    ent.debug_print_base( os );
    return( os );
}


void MyDXFEntity::process_group( class MyDXFFile *dxf )
{
    if( dxf->group_get_code() == 5 ) {
	_handle = dxf->group_get_string();
    } else if( dxf->group_get_code() == 8 ) {
	_layer = dxf->group_get_string();
    }
}




/* ************************************************************************** *
 * DXFLine                                                                    *
 * ************************************************************************** */


MyDXFLine::MyDXFLine( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity LINE\n";
#endif

    // Default values
    _p1[0] = _p1[1] = _p1[2] = 0.0;
    _p2[0] = _p2[1] = _p2[2] = 0.0;

    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 10 )
	    _p1[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 )
	    _p1[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 )
	    _p1[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 11 )
	    _p2[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 21 )
	    _p2[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 31 )
	    _p2[2] = dxf->group_get_double();
	else
	    process_group( dxf );
    }

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFLine::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		      const Transformation *t, const double range[4] ) const
{
    Vec4D x1 = t->transform( _p1 );
    Vec4D x2 = t->transform( _p2 );

    cairo_move_to( cairo, x1[0], x1[1] );
    cairo_line_to( cairo, x2[0], x2[1] );

    cairo_stroke( cairo );
}


void MyDXFLine::get_bbox( Vec3D &min, Vec3D &max, 
			  const class MyDXFFile *dxf, const Transformation *t ) const
{
    Vec4D x1 = t->transform( _p1 );
    Vec4D x2 = t->transform( _p2 );

    for( int a = 0; a < 3; a++ ) {
	if( x1[a] < x2[a] ) {
	    min[a] = x1[a];
	    max[a] = x2[a];
	} else {
	    min[a] = x2[a];
	    max[a] = x1[a];
	}
    }	

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "Line bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFLine::scale( class MyDXFFile *dxf, double s  )
{
    _p1 *= s;
    _p2 *= s;
}


int MyDXFLine::ray_cross( double x, double y ) const
{
    if( (x > _p1[0] && x < _p2[0]) || 
	(x < _p1[0] && x > _p2[0]) ) {

	// Calculate crossing y-coordinate.
	double t = (x-_p1[0])/(_p2[0]-_p1[0]);
        double cy = (1.0-t)*_p1[1] + t*_p2[1];
        // Boundary case y == cy is considered crossing.
        if( y >= cy )
            return( 1 );
    }

    // Exact crossing.
    if( x == _p1[0] && y >= _p1[1] )
        return( 2 );
    
    // No crossing.
    return( 0 );
}


void MyDXFLine::debug_print( std::ostream &os ) const
{
    std::cout << "LINE\n";
    std::cout << "  p1 = {" 
	      << _p1[0] << ", " 
	      << _p1[1] << ", " 
	      << _p1[2] << "}\n";
    std::cout << "  p2 = {" 
	      << _p2[0] << ", " 
	      << _p2[1] << ", " 
	      << _p2[2] << "}\n";
}


void MyDXFLine::set_start( const Vec3D &s )
{
    _p1 = s;
}


void MyDXFLine::set_end( const Vec3D &e )
{
    _p2 = e;
}


bool MyDXFLine::geom_same( const MyDXFLine &line, double eps ) const
{
    return( (norm2( _p1 - line._p1 ) < eps &&
	     norm2( _p2 - line._p2 ) < eps) ||
	    (norm2( _p2 - line._p1 ) < eps &&
	     norm2( _p1 - line._p2 ) < eps) );
}




/* ************************************************************************** *
 * DXFLWPolyline                                                              *
 * ************************************************************************** */


MyDXFLWPolyline::MyDXFLWPolyline( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity LWPOLYLINE\n";
#endif

    int32_t N = 0;
    int32_t i = -1;
    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 10 ) {
	    i++;
	    if( i == N )
		throw Error( ERROR_LOCATION, "Incorrect amount of vertices in LWPOLYLINE"
			     " on line " + dxf->linec() );
	    _p[i][0] = dxf->group_get_double();
	} else if( dxf->group_get_code() == 20 )
	    _p[i][1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 42 )
	    _p[i][2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 70 )
	    _flags = dxf->group_get_int16();
	else if( dxf->group_get_code() == 90 ) {
	    // Resize to number of vertices.
	    if( i != -1 )
		throw Error( ERROR_LOCATION, "Number of vertices specified after first "
			     "vertedin LWPOLYLINE on line " + dxf->linec() );
	    N = dxf->group_get_int32();
	    _p.resize( N );
	}

	else
	    process_group( dxf );
    }

    if( i != N-1 )
	throw Error( ERROR_LOCATION, "Incorrect amount of vertices in LWPOLYLINE"
		     " on line " + dxf->linec() );

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFLWPolyline::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
			    const Transformation *t, const double range[4] ) const
{
    if( _p.size() == 0 )
	return;

    Vec4D x = t->transform( _p[0] );
    cairo_move_to( cairo, x[0], x[1] );
    for( uint32_t a = 1; a < _p.size(); a++ ) {
	x = t->transform( _p[a] );
	cairo_line_to( cairo, x[0], x[1] );
    }
    if( closed() ) {
	x = t->transform( _p[0] );
	cairo_line_to( cairo, x[0], x[1] );
    }

    cairo_stroke( cairo );
}


void MyDXFLWPolyline::get_bbox( Vec3D &min, Vec3D &max, 
				const class MyDXFFile *dxf, const Transformation *t ) const
{
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    for( uint32_t a = 0; a < _p.size(); a++ ) {
	Vec4D x = t->transform( _p[a] );
	bbox_ppoint( min, max, x );
    }

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "Polyline bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFLWPolyline::scale( class MyDXFFile *dxf, double s  )
{
    for( uint32_t a = 0; a < _p.size(); a++ )
	_p[a] *= s;
}


int MyDXFLWPolyline::ray_cross( double x, double y ) const
{
#ifdef MYDXF_DEBUG
    std::cout << "ray_cross( x = " << x << ", y = " << y << " )\n";
    for( uint32_t a = 0; a < _p.size(); a++ )
	std::cout << "  p[" << a << "] = " << _p[a] << "\n";
    std::cout << "\n";
#endif

    Vec3D p1 = _p[0];

    // Go through all lines 
    int c = 0;
    for( uint32_t a = 1; a < _p.size(); a++ ) {
	Vec3D p2 = _p[a];

#ifdef MYDXF_DEBUG
	std::cout << "  Test against:\n";
	std::cout << "  p1 = " << p1 << "\n";
	std::cout << "  p2 = " << p2 << "\n";
#endif

	if( x == p1[0] && y >= p1[1] ) {
	    // Exact crossing.
#ifdef MYDXF_DEBUG
	    std::cout << "  Exact crossing\n";
#endif
	    return( 2 );
	} else if( (x > p1[0] && x < p2[0]) || 
		   (x < p1[0] && x > p2[0]) ) {

	    // Calculate crossing y-coordinate.
	    double t = (x-p1[0])/(p2[0]-p1[0]);
	    double cy = (1.0-t)*p1[1] + t*p2[1];
	    // Boundary case y == cy is considered crossing.
	    if( y >= cy ) {
#ifdef MYDXF_DEBUG
		std::cout << "  Crossing\n";
#endif
		c = !c;
	    } else {
#ifdef MYDXF_DEBUG
		std::cout << "  No crossing\n";
#endif
	    }
	    
	}

	p1 = p2;
    }

    return( c );
}


void MyDXFLWPolyline::debug_print( std::ostream &os ) const
{
    std::cout << "LWPOLYLINE\n";
    std::cout << "  flags = " << _flags << "\n"; 
    for( uint32_t a = 0; a < _p.size(); a++ ) {
	std::cout << "  p[" << a << "] = {" 
		  << _p[a][0] << ", " 
		  << _p[a][1] << ", " 
		  << _p[a][2] << "}\n";
    }
}


void MyDXFLWPolyline::set_start( const Vec3D &s )
{
    _p[0] = s;
}


void MyDXFLWPolyline::set_end( const Vec3D &e )
{
    _p[_p.size()-1] = e;
}


bool MyDXFLWPolyline::geom_same( const MyDXFLWPolyline &line, double eps ) const
{
    if( _p.size() != line._p.size() )
	return( 1 );

    for( uint32_t a = 0; a < _p.size(); a++ ) {
	if( norm2( _p[a] - line._p[a] ) < eps )
	    return( 1 );
    }

    return( 0 );
}




/* ************************************************************************** *
 * DXFArc                                                                     *
 * ************************************************************************** */


MyDXFCircle::MyDXFCircle( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity CIRCLE\n";
#endif

    // Default values
    _pc[0] = _pc[1] = _pc[2] = 0.0;
    _r = 1.0;

    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 10 )
	    _pc[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 )
	    _pc[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 )
	    _pc[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 40 )
	    _r = dxf->group_get_double();

	else
	    process_group( dxf );
    }

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFCircle::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
			const Transformation *t, const double range[4] ) const
{
    std::vector<double> cx;
    std::vector<double> cy;

    // Divide circle to linear pieces and keep length of linear pieces
    // to less than 4 px

    // Start with three points
    double cf = 2.0*M_PI/3.0;
    for( int i = 0; i < 3; i++ ) {
	double a = i*cf;
	Vec4D p = t->transform( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );
	cx.push_back( p[0] );
	cy.push_back( p[1] );
    }

    // Loop until error is small enough
    for( int iter = 0; iter < MAX_ARC_CHOP; iter++ ) {

	//std::cout << "LOOP: N = " << cx.size() << "\n";
	//for( int i = 0; i < (int)cx.size(); i++ )
	//std::cout << "  " << cx[i] << "\t" << cy[i] << "\n";

	// Calculate maximum step length
	double mstep = 0.0;
	for( int i = 0; i < (int)cx.size()-1; i++ ) {
	    double sx = cx[i+1] - cx[i];
	    double sy = cy[i+1] - cy[i];
	    double step = sx*sx + sy*sy;
	    if( step > mstep )
		mstep = step;
	}
	// last gap
	double sx = cx[0] - cx[cx.size()-1];
	double sy = cy[0] - cy[cx.size()-1];
	double step = sx*sx + sy*sy;
	if( step > mstep )
	    mstep = step;

	// Calculate maximum error between most recent points and the
	// linear interpolation from the last round
	/*
	double maxerr = 0.0;
	for( int i = 1; i < (int)cx.size(); i += 2 ) {
	    double tx = 0.5*(cx[i-1] + cx[i+1]);
	    double ty = 0.5*(cy[i-1] + cy[i+1]);
	    tx -= cx[i];
	    ty -= cy[i];
	    double err = tx*tx + ty*ty; // square of error
	    if( err > maxerr )
		maxerr = err;
	}
	*/

	//std::cout << "mstep = " << mstep << "\n";

	// If step < 4.0
	if( mstep < 16.0 )
	    break;

	// Add more points
	cx.resize( 2*cx.size() );
	cy.resize( 2*cy.size() );
	double cf = 2.0*M_PI/cx.size();
	for( int i = cx.size()-1; i > 0; i-- ) {
	    
	    // New point
	    double a = i*cf;
	    Vec4D p = t->transform( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );
	    cx[i] = p[0];
	    cy[i] = p[1];
	    i--;

	    // Old point
	    cx[i] = cx[i/2];
	    cy[i] = cy[i/2];
	}

	//std::cout << "\n";
    }

    // Error small enough -> draw lines
    cairo_move_to( cairo, cx[0], cy[0] );
    for( int i = 0; i < (int)cx.size(); i++ )
	cairo_line_to( cairo, cx[i], cy[i] );
    cairo_line_to( cairo, cx[0], cy[0] );
    cairo_stroke( cairo );
}


void MyDXFCircle::get_bbox( Vec3D &min, Vec3D &max, 
			    const class MyDXFFile *dxf, const Transformation *t ) const
{
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    // Chop circle to 8 pieces for bbox
    double cf = 2.0*M_PI/8.0;
    for( int i = 0; i < 8; i++ ) {
	double a = i*cf;
	Vec3D p( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] );
	Vec4D x = t->transform( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );
	bbox_ppoint( min, max, x );
    }

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "Circle bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFCircle::scale( class MyDXFFile *dxf, double s )
{
    _pc *= s;
    _r  *= s;
}


int MyDXFCircle::ray_cross( double x, double y ) const
{
    if( (x <= _pc[0]-_r || x >= _pc[0]+_r) )
	// No crossing
	return( 0 );

    double b = x-_pc[0];
    double t = sqrt( _r*_r - b*b );
    double cy1 = _pc[1] + t;
    double cy2 = _pc[1] - t;

    if( y <= cy2 )
	return( 0 );
    else if( y >= cy1 )
	return( 0 );

    return( 1 );
}


void MyDXFCircle::debug_print( std::ostream &os ) const
{
    std::cout << "CIRCLE\n";
    std::cout << "  p = {" 
	      << _pc[0] << ", " 
	      << _pc[1] << ", " 
	      << _pc[2] << "}\n";
    std::cout << "  r = " << _r << "\n";
}


bool MyDXFCircle::geom_same( const MyDXFCircle &circle, double eps ) const
{
    return( norm2( _pc - circle._pc ) < eps &&
	    fabs( _r - circle._r ) < eps );
}




/* ************************************************************************** *
 * DXFArc                                                                     *
 * ************************************************************************** */


MyDXFArc::MyDXFArc( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity ARC\n";
#endif

    // Default values
    _pc[0] = _pc[1] = _pc[2] = 0.0;
    _r = 1.0;
    _ang1 = 0.0;
    _ang2 = 2.0*M_PI;

    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 10 )
	    _pc[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 )
	    _pc[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 )
	    _pc[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 40 )
	    _r = dxf->group_get_double();

	else if( dxf->group_get_code() == 50 ) {
	    double ang1 = M_PI*dxf->group_get_double()/180.0;
	    // Enforce between 0 and 2 pi
	    _ang1 = ang1 - 2.0*M_PI*floor( ang1/(2.0*M_PI) );
	} else if( dxf->group_get_code() == 51 ) {
	    double ang2 = M_PI*dxf->group_get_double()/180.0;
	    // Enforce between 0 and 2 pi
	    _ang2 = ang2 - 2.0*M_PI*floor( ang2/(2.0*M_PI) );
	}

	else
	    process_group( dxf );
    }

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFArc::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		     const Transformation *t, const double range[4] ) const
{
    std::vector<double> cx;
    std::vector<double> cy;

    // Divide arc to linear pieces and keep length of linear pieces
    // to less than 4 px

    //std::cout << "\n\nARC plot\n";
    //std::cout << "  ang1 = " << _ang1 << "\n";
    //std::cout << "  ang2 = " << _ang2 << "\n";

    // Start with three points: ends and one in between
    double adiff = _ang2-_ang1;
    if( adiff < 0.0 )
	adiff = 2.0*M_PI + adiff;
    //std::cout << "  adiff = " << adiff << "\n";

    for( int i = 0; i < 3; i++ ) {
	double a = _ang1 + i*adiff/2.0;
	if( a > 2.0*M_PI )
	    a -= 2.0*M_PI;
	Vec4D p = t->transform( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );
	//std::cout << "  a = " << a << "\n";
	//std::cout << "  p = " << p[0] << "\t" << p[1] << "\n";
	cx.push_back( p[0] );
	cy.push_back( p[1] );
    }

    // Loop until error is small enough
    for( int iter = 0; iter < MAX_ARC_CHOP; iter++ ) {

	//std::cout << "LOOP: N = " << cx.size() << "\n";
	//for( int i = 0; i < (int)cx.size(); i++ )
	//    std::cout << "  " << cx[i] << "\t" << cy[i] << "\n";

	// Calculate maximum step length
	double mstep = 0.0;
	for( int i = 0; i < (int)cx.size()-1; i++ ) {
	    double sx = cx[i+1] - cx[i];
	    double sy = cy[i+1] - cy[i];
	    double step = sx*sx + sy*sy;
	    if( step > mstep )
		mstep = step;
	}
	//std::cout << "mstep = " << mstep << "\n";
	//std::cout << "\n";

	// If step < 4.0
	if( mstep < 16.0 )
	    break;

	// Add more points
	cx.resize( 2*cx.size()-1 );
	cy.resize( 2*cy.size()-1 );
	double cf = adiff/(cx.size()-1);
	for( int i = cx.size()-1; i > 1; i-- ) {

	    // Old point
	    cx[i] = cx[i/2];
	    cy[i] = cy[i/2];
	    i--;

	    // New point
	    double a = _ang1 + i*cf;
	    if( a > 2.0*M_PI )
		a -= 2.0*M_PI;
	    Vec4D p = t->transform( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );
	    //std::cout << "  a = " << a << "\n";
	    //std::cout << "  p = " << p[0] << "\t" << p[1] << "\n";
	    cx[i] = p[0];
	    cy[i] = p[1];
	}
    }

    // Error small enough -> draw lines
    cairo_move_to( cairo, cx[0], cy[0] );
    for( int i = 0; i < (int)cx.size(); i++ )
	cairo_line_to( cairo, cx[i], cy[i] );
    cairo_stroke( cairo );
}


void MyDXFArc::get_bbox( Vec3D &min, Vec3D &max, 
			 const class MyDXFFile *dxf, const Transformation *t ) const
{
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    // Chop arc to 8 pieces for bbox
    double adiff = _ang2-_ang1;
    if( adiff < 0.0 )
	adiff = 2.0*M_PI + adiff;

    for( int i = 0; i < 8; i++ ) {
	double a = _ang1 + i*adiff/7.0;
	if( a > 2.0*M_PI )
	    a -= 2.0*M_PI;
	Vec4D x = t->transform( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );
	bbox_ppoint( min, max, x );
    }

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "Arc bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFArc::scale( class MyDXFFile *dxf, double s )
{
    _pc *= s;
    _r  *= s;
}


int MyDXFArc::ray_cross( double x, double y ) const
{
    if( (x <= _pc[0]-_r || x >= _pc[0]+_r) )
	// No crossing
	return( 0 );

    double dx = x-_pc[0];
    double dy = sqrt( _r*_r - dx*dx );
    double cy1 = _pc[1] + dy;
    double cy2 = _pc[1] - dy;

    if( y <= cy2 )
	// Under the circle -> no crossing
	return( 0 );

    double a1 = atan2( dy,dx );
    if( a1 < 0.0 )
	a1 += 2.0*M_PI;

    double a2 = atan2( -dy,dx );
    if( a2 < 0.0 )
	a2 += 2.0*M_PI;

    int c = 0;
    if( _ang1 < _ang2 ) {
	if( a1 > _ang1 && a1 < _ang2 && y > cy1 )
	    c = !c;
	if( a2 > _ang1 && a2 < _ang2 && y > cy2 )
	    c = !c;
    } else {
	if( (a1 < _ang2 || a1 > _ang1) && y > cy1 )
	    c = !c;
	if( (a2 < _ang2 || a2 > _ang1) && y > cy2 )
	    c = !c;
    }

    return( c );
}

void MyDXFArc::set_ang1( double ang1 ) 
{
    // Enforce between 0 and 2 pi
    _ang1 = ang1 - 2.0*M_PI*floor( ang1/(2.0*M_PI) );
}


void MyDXFArc::set_ang2( double ang2 )
{ 
    // Enforce between 0 and 2 pi
    _ang2 = ang2 - 2.0*M_PI*floor( ang2/(2.0*M_PI) );
}


void MyDXFArc::debug_print( std::ostream &os ) const
{
    std::cout << "ARC\n";
    std::cout << "  p = {" 
	      << _pc[0] << ", " 
	      << _pc[1] << ", " 
	      << _pc[2] << "}\n";
    std::cout << "  r = " << _r << "\n";
    std::cout << "  ang1 = " << _ang1 << "\n";
    std::cout << "  ang2 = " << _ang2 << "\n";
}


bool MyDXFArc::geom_same( const MyDXFArc &arc, double eps ) const
{
    return( norm2( _pc - arc._pc ) < eps &&
	    fabs( _r - arc._r ) < eps &&
	    norm2( start() - arc.start() ) < eps &&
	    norm2( end() - arc.end() ) < eps );
}


void MyDXFArc::set_start( const Vec3D &s )
{
    Vec3D e = end();
    set_center_point( s, e );
}


void MyDXFArc::set_end( const Vec3D &e )
{
    Vec3D s = start();
    set_center_point( s, e );
}


void MyDXFArc::set_center_point( const Vec3D &s, const Vec3D &e )
{
    if( norm2( e-s ) < 2.0*_r ) {
	Vec3D u = 0.5*(e-s);
	Vec3D v( -u[1], u[0] );

	double A = v[0]*v[0] + v[1]*v[1];
	double B = 2.0*( v[0]*u[0] + v[1]*u[1] );
	double C = u[0]*u[0] + u[1]*u[1] - _r*_r;

	double t1, t2;
	uint32_t n = solve_quadratic( A, B, C, &t1, &t2 );
	if( n != 2 )
	    throw Error( ERROR_LOCATION, "Less than two roots found" );
	
	// solver gives solutions in ascending order and t2 = -t1
	_pc = s+u+t2*v;
    } else {
	_pc = 0.5*(s+e);
	_r = 0.5*norm2( e-s );
    }

    Vec3D dd = s-_pc;
    double dx = dd[0];
    double dy = dd[1];
    _ang1 = atan2( dy,dx );
    if( _ang1 < 0.0 )
	_ang1 += 2.0*M_PI;

    dd = e-_pc;
    dx = dd[0];
    dy = dd[1];
    _ang2 = atan2( dy,dx );
    if( _ang2 < 0.0 )
	_ang2 += 2.0*M_PI;
}




/* ************************************************************************** *
 * DXFMText                                                                   *
 * ************************************************************************** */


MyDXFMText::MyDXFMText( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity MTEXT\n";
#endif

    // Default values
    _p[0] = _p[1] = _p[2] = 0.0;
    _text_height = 1.0;
    _rect_width = 1.0;
    _attachment_point = ATTACHMENT_POINT_TOP_LEFT;
    _drawing_direction = DRAWING_DIRECTION_LEFT_TO_RIGHT;

    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 10 )
	    _p[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 )
	    _p[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 )
	    _p[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 40 )
	    _text_height = dxf->group_get_double();
	else if( dxf->group_get_code() == 41 )
	    _rect_width = dxf->group_get_double();

	else if( dxf->group_get_code() == 71 )
	    _attachment_point = dxf->group_get_int16();
	else if( dxf->group_get_code() == 72 )
	    _drawing_direction = dxf->group_get_int16();

	else if( dxf->group_get_code() == 3 || dxf->group_get_code() == 1 )
	    _text += dxf->group_get_string();
	
	else
	    process_group( dxf );
    }

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFMText::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		       const Transformation *t, const double range[4] ) const
{

}


void MyDXFMText::get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const
{
    std::cout << "Warning: bounding box for MText entity not implemented\n";
    min = _p;
    max = _p;

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "MText bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFMText::scale( class MyDXFFile *dxf, double s )
{
    _p *= s;
    _text_height *= s;
    _rect_width *= s;
}


void MyDXFMText::debug_print( std::ostream &os ) const
{
    std::cout << "MTEXT\n";
    std::cout << "  p = {" 
	      << _p[0] << ", " 
	      << _p[1] << ", " 
	      << _p[2] << "}\n";
    std::cout << "  text_height = " << _text_height << "\n"; 
    std::cout << "  rect_width = " << _rect_width << "\n"; 
    std::cout << "  attachment_point = " << _attachment_point << "\n"; 
    std::cout << "  drawing_direction = " << _drawing_direction << "\n"; 
}




/* ************************************************************************** *
 * DXFMInsert                                                                 *
 * ************************************************************************** */


MyDXFInsert::MyDXFInsert()
    : _rotation(0.0), _col_count(1), _row_count(1), _col_spacing(0.0), _row_spacing(0.0)
{
    // Default values
    _scale[0] = _scale[1] = _scale[2] = 1.0;
}


MyDXFInsert::MyDXFInsert( class MyDXFFile *dxf )
    : _rotation(0.0), _col_count(1), _row_count(1), _col_spacing(0.0), _row_spacing(0.0)
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity INSERT\n";
#endif

    // Default values
    _scale[0] = _scale[1] = _scale[2] = 1.0;

    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 2 )
	    _block_name = dxf->group_get_string();

	else if( dxf->group_get_code() == 10 )
	    _p[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 )
	    _p[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 )
	    _p[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 41 )
	    _scale[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 42 )
	    _scale[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 43 )
	    _scale[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 44 )
	    _col_spacing = dxf->group_get_double();
	else if( dxf->group_get_code() == 45 )
	    _row_spacing = dxf->group_get_double();

	else if( dxf->group_get_code() == 50 ) {
	    double rot = M_PI*dxf->group_get_double()/180.0;
	    // Enforce between 0 and 2 pi
	    _rotation = rot - 2.0*M_PI*floor( rot/(2.0*M_PI) );
	}

	else if( dxf->group_get_code() == 70 )
	    _col_count = dxf->group_get_int16();
	else if( dxf->group_get_code() == 71 )
	    _row_count = dxf->group_get_int16();

	else
	    process_group( dxf );
    }

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFInsert::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
			const Transformation *t, const double range[4] ) const
{
    // Fetch block data
    const MyDXFBlocks *blocks = dxf->get_blocks();
    const MyDXFBlock *b = blocks->get_by_name( _block_name );
    if( !b )
	return;

    Transformation t2 = *t;
    t2.translate( _p );
    t2.scale( _scale );
    t2.rotate_z( _rotation );

    for( int16_t col = 0; col < _col_count; col++ ) {
	for( int16_t row = 0; row < _row_count; row++ ) {
	    Transformation t3 = t2;
	    t3.translate( Vec3D( col*_col_spacing, row*_row_spacing, 0.0 ) );
	    b->plot( dxf, cairo, &t3, range );
	}
    }
}


void MyDXFInsert::get_bbox( Vec3D &min, Vec3D &max, 
			    const class MyDXFFile *dxf, const Transformation *t ) const
{
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    // Fetch block data
    const MyDXFBlocks *blocks = dxf->get_blocks();
    const MyDXFBlock *b = blocks->get_by_name( _block_name );
    if( !b )
	return;

    Transformation t2 = *t;
    t2.translate( _p );
    t2.scale( _scale );
    t2.rotate_z( _rotation );

    for( int16_t col = 0; col < _col_count; col++ ) {
	for( int16_t row = 0; row < _row_count; row++ ) {
	    Transformation t3 = t2;
	    t3.translate( Vec3D( col*_col_spacing, row*_row_spacing, 0.0 ) );
	    Vec3D mi, ma;
	    b->get_bbox( mi, ma, dxf, &t3 );

	    for( int b = 0; b < 3; b++ ) {
		if( mi[b] < min[b] )
		    min[b] = mi[b];
		if( ma[b] > max[b] )
		    max[b] = ma[b];
	    }
	}
    }

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "Insert bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFInsert::scale( class MyDXFFile *dxf, double s )
{
    // Fetch block data
    MyDXFBlocks *blocks = dxf->get_blocks();
    MyDXFBlock *b = blocks->get_by_name( _block_name );
    if( !b )
	return;

    b->scale( dxf, s );
}


void MyDXFInsert::debug_print( std::ostream &os ) const
{
    std::cout << "INSERT\n";
    std::cout << "  block_name = \'" << _block_name << "\'\n";
    std::cout << "  p = {" 
	      << _p[0] << ", " 
	      << _p[1] << ", " 
	      << _p[2] << "}\n";
    std::cout << "  scale = {" 
	      << _scale[0] << ", " 
	      << _scale[1] << ", " 
	      << _scale[2] << "}\n";
    std::cout << "  rotation = " << _rotation << "\n";
    std::cout << "  col_count = " << _col_count << "\n";
    std::cout << "  row_count = " << _row_count << "\n";
    std::cout << "  col_spacing = " << _col_spacing << "\n";
    std::cout << "  row_spacing = " << _row_spacing << "\n";
}




/* ************************************************************************** *
 * DXFEntities                                                                *
 * ************************************************************************** */


MyDXFEntities::MyDXFEntities( MyDXFEntities *ent, MyDXFEntitySelection *sel )
{
    for( size_t a = 0; a < sel->size(); a++ ) {
	MyDXFEntity *e = ent->get_entity( (*sel)(a) );
	_entities.push_back( e->copy() );
    }
}


MyDXFEntities::MyDXFEntities( class MyDXFFile *dxf, bool reading_blocks )
{
#ifdef MYDXF_DEBUG
    std::cout << "Reading ENTITIES\n";
#endif

    // Title of first entity already read if reading BLOCKS section
    if( !reading_blocks )
	dxf->read_group();

    while( dxf->group_get_code() != -1 ) {

	if( dxf->group_get_code() != 0 ) {
	    dxf->read_group();
	    continue; // Skip unknown input
	}
	else if( dxf->group_get_string() == "ENDSEC" ||
		 dxf->group_get_string() == "ENDBLK" ) {
	    break; // Done with entities
	}

	// Read entity type
	if( dxf->group_get_string() == "LINE" ) {
	    _entities.push_back( new MyDXFLine( dxf ) );
	} else if( dxf->group_get_string() == "CIRCLE" ) {
	    _entities.push_back( new MyDXFCircle( dxf ) );
	} else if( dxf->group_get_string() == "MTEXT" ) {
	    _entities.push_back( new MyDXFMText( dxf ) );
	} else if( dxf->group_get_string() == "INSERT" ) {
	    _entities.push_back( new MyDXFInsert( dxf ) );
	} else if( dxf->group_get_string() == "ARC" ) {
	    _entities.push_back( new MyDXFArc( dxf ) );
	} else if( dxf->group_get_string() == "LWPOLYLINE" ) {
	    _entities.push_back( new MyDXFLWPolyline( dxf ) );
	} else {
	    if( dxf->wlevel() )
		std::cout << "Skipping unknown entity \'" << dxf->group_get_string() << "\'\n";
	    dxf->read_group();
	}
    }
}


MyDXFEntities::~MyDXFEntities()
{
    // Free data
    for( size_t a = 0; a < _entities.size(); a++ )
	delete _entities[a];
}


void MyDXFEntities::write( class MyDXFFile *dxf, std::ofstream &_ostr )
{
    dxf->write_group( 0, "SECTION" );
    dxf->write_group( 2, "ENTITIES" );

    for( size_t a = 0; a < _entities.size(); a++ )
	

    dxf->write_group( 0, "ENDSEC" );
}


bool MyDXFEntities::inside_loop( MyDXFEntitySelection *selection, double x, double y, double eps )
{
    for( uint32_t b = 0; b < 2; b++ ) {

	int stat = 0;
        int par = 0;
        for( uint32_t a = 0; a < selection->size(); a++ ) {
	    MyDXFEntity *e = _entities[(*selection)(a)];
	    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	    if( !pe )
		throw Error( ERROR_LOCATION, "Not a path entity" );
	    stat = pe->ray_cross( x, y );
	    if( stat == 1 ) 
		par = !par;
	    else if( stat == 2 ) 
		break;
        }
        
        if( stat != 2 )
            return( par );

	// Perturbation
        x += eps;
    }

    throw Error( ERROR_LOCATION, "Perturbation failed" );
}


MyDXFEntitySelection *MyDXFEntities::selection_layer( const std::string &layername ) const
{
    MyDXFEntitySelection *selection = new MyDXFEntitySelection();

    for( size_t a = 0; a < _entities.size(); a++ ) {
	if( _entities[a]->get_layer() == layername )
	    selection->add_entity( a );
    }

    return( selection );
}


MyDXFEntitySelection *MyDXFEntities::selection_all( void ) const
{
    MyDXFEntitySelection *selection = new MyDXFEntitySelection();

    for( size_t a = 0; a < _entities.size(); a++ ) {
	selection->add_entity( a );
    }

    return( selection );
}


std::ostream &operator<<( std::ostream &os, const MyDXFEntitySelection &sel )
{
    if( sel._selection.size() > 1 ) {
	os << "{";
	for( size_t a = 0; a < sel._selection.size()-1; a++ )
	    os << sel._selection[a] << ", ";
	os << sel._selection[sel._selection.size()-1] << "}";
    } else {
	os << "{" << sel._selection[0] << "}";
    }
    
    return( os );
}


bool MyDXFEntities::geom_same( uint32_t a, uint32_t b, double eps ) const
{
    // Test lines
    const MyDXFLine *line1 = dynamic_cast<MyDXFLine *>( _entities[a] );
    const MyDXFLine *line2 = dynamic_cast<MyDXFLine *>( _entities[b] );
    if( line1 && line2 )
	return( line1->geom_same( *line2, eps ) );

    // Test arcs
    const MyDXFArc *arc1 = dynamic_cast<MyDXFArc *>( _entities[a] );
    const MyDXFArc *arc2 = dynamic_cast<MyDXFArc *>( _entities[b] );
    if( arc1 && arc2 )
	return( arc1->geom_same( *arc2, eps ) );
    
    // Test circles
    const MyDXFCircle *circle1 = dynamic_cast<MyDXFCircle *>( _entities[a] );
    const MyDXFCircle *circle2 = dynamic_cast<MyDXFCircle *>( _entities[b] );
    if( circle1 && circle2 )
	return( circle1->geom_same( *circle2, eps ) );

    // Test lwpolylines
    const MyDXFLWPolyline *lwpline1 = dynamic_cast<MyDXFLWPolyline *>( _entities[a] );
    const MyDXFLWPolyline *lwpline2 = dynamic_cast<MyDXFLWPolyline *>( _entities[b] );
    if( lwpline1 && lwpline2 )
	return( lwpline1->geom_same( *lwpline2, eps ) );

    return( false );
}


MyDXFEntitySelection *MyDXFEntities::selection_path_loop( MyDXFEntitySelection *selection,
							  double eps )
{
    uint32_t a;
    MyDXFEntitySelection *subsel = new MyDXFEntitySelection();

    std::vector<bool>     stdir;     // Direction of entities on stack
    std::vector<uint32_t> stack;     // Stack of entities for loop
    bool done[selection->size()];    // List of processed entities

    // Initialize
    for( a = 0; a < selection->size(); a++ ) {
	done[a] = false;
    }

    /*
    std::cout << "*************************************************************\n";
    std::cout << "selection = {";
    for( a = 0; a < selection->size()-1; a++ )
	std::cout << (*selection)(a) << ", ";
    std::cout << (*selection)(a) << "}";
    */

    // Remove non-path objects and process self-looped path objects
    for( a = 0; a < selection->size(); a++ ) {
	MyDXFEntity *e = _entities[(*selection)(a)];
	MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	if( !pe )
	    done[a] = true;
	else if( pe->start() == pe->end() ) {
	    // Add entity to list and remove duplicates
	    //std::cout << "Adding " << (*selection)(a) << " to selection";
	    subsel->add_entity( (*selection)(a) );
	    done[a] = true;
	    for( uint32_t b = a+1; b < selection->size(); b++ ) {
		MyDXFEntity *eb = _entities[(*selection)(b)];
		MyDXFPathEntity *peb = dynamic_cast<MyDXFPathEntity *>( eb );
		if( peb && peb->start() == peb->end() ) {
		    if( geom_same( (*selection)(a), (*selection)(b), eps ) ) {
			done[b] = true;
			continue;
		    }
		}
	    }
	}
    }

    // Loop until all entities done
    while( 1 ) {

#ifdef MYDXF_DEBUG
	std::cout << "\n\ndone = {";
	for( a = 0; a < selection->size()-1; a++ )
	    std::cout << done[a] << ", ";
	std::cout << done[a] << "}\n";
#endif

	if( stack.size() == 0 ) {
	    // Start with first unprocessed of the selected entities,
	    // which is a path object
	    for( a = 0; a < selection->size(); a++ ) {
		if( !done[a] )
		    break;
	    }
	    if( a == selection->size() )
		break; // No entities left
	    // Add to stack
	    done[a] = true;
	    stdir.push_back( true );
	    stack.push_back( (*selection)(a) );

#ifdef MYDXF_DEBUG
	    std::cout << "Starting stack with " << (*selection)(a) << "\n";
#endif
	}

	// Check if loop done, check if endpoint of last entity
	// on stack matches some starting point on stack
	Vec3D end, start;
	if( stdir.back() ) {
	    MyDXFEntity *e = _entities[ stack.back() ];
	    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	    end = pe->end();
	} else {
	    MyDXFEntity *e = _entities[ stack.back() ];
	    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	    end = pe->start();
	}
#ifdef MYDXF_DEBUG
	std::cout << "Check if we have loop m. "
		  << "end = " << end << "\n";
#endif

	for( a = 0; a < stack.size(); a++ ) {
	    if( stdir[a] ) {
		MyDXFEntity *e = _entities[ stack[a] ];
		MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
		start = pe->start();
	    } else {
		MyDXFEntity *e = _entities[ stack[a] ];
		MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
		start = pe->end();
	    }
#ifdef MYDXF_DEBUG
	    std::cout << "  Stack[" << a << "]: "
		      << "stdir = " << stdir[a] << ", "
		      << "start = " << start << "\n";
#endif
	    // If loop closed
	    if( norm2( start-end ) < eps ) {

#ifdef MYDXF_DEBUG
		std::cout << "  Match found\n";
#endif
		// Fix all ends from stack and add to subsel list
		for( int32_t b = stack.size()-1; b >= (int32_t)a; b-- ) {

		    uint32_t c = stack.back();
		    MyDXFEntity *e = _entities[c];
		    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );

		    if( stdir.back() )
			pe->set_end( start );
		    else
			pe->set_start( start );

		    if( stdir[b] )
			start = pe->start();
		    else
			start = pe->end();
		    
		    subsel->add_entity( c );
		    stdir.pop_back();
		    stack.pop_back();
		}
		break;
	    }
	}
	if( stack.size() == 0 )
	    continue;

	// Search for an unprocessed and selected entity, which is a
	// path object and start matches with the end of the last
	// entity on stack
	if( stdir.back() ) {
	    MyDXFEntity *e = _entities[ stack.back() ];
	    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	    end = pe->end();
	} else {
	    MyDXFEntity *e = _entities[ stack.back() ];
	    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	    end = pe->start();
	}
#ifdef MYDXF_DEBUG
	std::cout << "Search for ent. matching end = "
		  << end << "\n";
#endif

	for( a = 0; a < selection->size(); a++ ) {

	    if( done[a] || (*selection)(a) == stack.back() ) 
		continue;

	    MyDXFEntity *e = _entities[ (*selection)(a) ];
	    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	    start = pe->start();

#ifdef MYDXF_DEBUG
	    std::cout << "  Entity " << std::setw(4) << (*selection)(a) << ": "
		      << "        start = " << start << "\n";
#endif
	    if( norm2( start-end ) < eps ) {
#ifdef MYDXF_DEBUG
		std::cout << "  Match\n";
#endif
		// Check if entities are same within eps
		if( geom_same( stack.back(), (*selection)(a), eps ) ) {
#ifdef MYDXF_DEBUG
		    std::cout << "  Removing duplicate\n";
#endif
		    // Remove duplicate
		    done[a] = true;
		    continue;
		}
		// Add to stack
		stdir.push_back( true );
		stack.push_back( (*selection)(a) );
		done[a] = true;
		break;
	    } else {

		start = pe->end();
#ifdef MYDXF_DEBUG
		std::cout << "  No match, other end: "
			  << "  end = " << start << "\n";
#endif

		if( norm2( start-end ) < eps ) {
#ifdef MYDXF_DEBUG
		    std::cout << "  Match\n";
#endif
		    // Check if entities are same within eps
		    if( geom_same( stack.back(), (*selection)(a), eps ) ) {
#ifdef MYDXF_DEBUG
			std::cout << "  Removing duplicate\n";
#endif
			// Remove duplicate
			done[a] = true;
			continue;
		    }
		    // Add to stack backwards
		    stdir.push_back( false );
		    stack.push_back( (*selection)(a) );
		    done[a] = true;
		    break;
		}
#ifdef MYDXF_DEBUG
		std::cout << "  No match\n";
#endif
	    }
	}
	if( a == selection->size() ) {
#ifdef MYDXF_DEBUG
	    std::cout << "  No matching entity found, removing last of stack\n";
#endif

	    // No matching entity found. Remove last entity of stack
	    // and mark it done
	    stdir.pop_back();
	    stack.pop_back();
	}

    }

    return( subsel );
}


void MyDXFEntities::plot( const MyDXFEntitySelection *selection, const class MyDXFFile *dxf, 
			  cairo_t *cairo, const Transformation *t, const double range[4] ) const
{
    if( selection ) {
	// Go through selection
	for( uint32_t a = 0; a < selection->size(); a++ ) {
	    MyDXFEntity *e = _entities[(*selection)(a)];
	    e->plot( dxf, cairo, t, range );
	}
    } else {
	// Plot all entities 
	for( uint32_t a = 0; a < _entities.size(); a++ ) {
	    MyDXFEntity *e = _entities[a];
	    e->plot( dxf, cairo, t, range );
	}
    }
}


void MyDXFEntities::get_bbox( const MyDXFEntitySelection *selection, Vec3D &min, Vec3D &max, 
			      const class MyDXFFile *dxf, const Transformation *t ) const
{
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    if( selection ) {
	for( size_t a = 0; a < selection->size(); a++ ) {
	    MyDXFEntity *e = _entities[(*selection)(a)];
	    Vec3D mi, ma;
	    e->get_bbox( mi, ma, dxf, t );

	    for( int b = 0; b < 3; b++ ) {
		if( mi[b] < min[b] )
		    min[b] = mi[b];
		if( ma[b] > max[b] )
		    max[b] = ma[b];
	    }
	}
    } else {
	for( size_t a = 0; a < _entities.size(); a++ ) {
	    MyDXFEntity *e = _entities[a];
	    Vec3D mi, ma;
	    e->get_bbox( mi, ma, dxf, t );
	    for( int b = 0; b < 3; b++ ) {
		if( mi[b] < min[b] )
		    min[b] = mi[b];
		if( ma[b] > max[b] )
		    max[b] = ma[b];
	    }
	}
    }

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "Entities bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFEntities::scale( MyDXFEntitySelection *selection, class MyDXFFile *dxf, double s )
{
    if( selection ) {
	// Go through selection
	for( size_t a = 0; a < selection->size(); a++ ) {
	    MyDXFEntity *e = _entities[(*selection)(a)];
	    e->scale( dxf, s );
	}    
    } else {
	// Scale all entities 
	for( size_t a = 0; a < _entities.size(); a++ ) {
	    MyDXFEntity *e = _entities[a];
	    e->scale( dxf, s );
	}
    }
}


void MyDXFEntities::debug_print( std::ostream &os ) const
{
    os << "*** Section ENTITIES **************************************\n";
    
    for( size_t a = 0; a < _entities.size(); a++ ) {
	MyDXFEntity *e = _entities[a];
	os << *e;
    }
    
    os << "\n";
}



