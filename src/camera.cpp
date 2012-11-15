/*! \file camera.cpp
 *  \brief Camera for 3D %Geometry viewer
 */

/* Copyright (c) 2012 Taneli Kalvas. All rights reserved.
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


#include "camera.hpp"


PerspectiveCamera::PerspectiveCamera()
    : _width(640), _height(480), _zoomFactor(1.0), 
      _near(1.0), _far(100.0), 
      _target(Vec3D(0,0,0)), _up(Vec3D(0,1,0)), _location(Vec3D(0,0,1))
{

}
    

PerspectiveCamera::~PerspectiveCamera()
{

}


void PerspectiveCamera::gl_initalize_camera( void )
{
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    if( _width > _height ) {
	_left = -1.0;
	_right = 1.0;
	_bottom = -(double)_height/_width;
	_top = (double)_height/_width;
    } else {
	_left = -(double)_width/_height;
	_right = (double)_width/_height;
	_bottom = -1.0;
	_top = 1.0;
    }
    glFrustum( _left*_zoom, 
               _right*_zoom, 
               _bottom*_zoom, 
               _top*_zoom, 
               _near, _far );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity( );
    gluLookAt( _camera[0],_camera[1],_camera[2], 
	       _target[0],_target[1],_target[2], 
	       _up[0],_up[1],_up[2] );
}


void PerspectiveCamera::set_size( double width, double height )
{
    _width = width;
    _height = height;
}


void PerspectiveCamera::set_view_relative( double x, double y, double fac )
{
    double u = _zoom*(_left + (_right-_left)*(x/_width));
    double v = _zoom*(_bottom + (_top-_bottom)*(1.0-y/_height));
    Vec3D viewdir = target-camera;
    viewdir.normalize();
    viewdir *= _near;
    Vec3D right = cross(viewdir,_up);
    Vec3D real_up = cross(right,viewdir);
    right.normalize();
    real_up.normalize();
    _target = _camera + viewdir + u*right + v*real_up;
    _zoom *= fac;
}


void PerspectiveCamera::set_target_location( const Vec3D &target )
{
    _target = target;
}


void PerspectiveCamera::set_camera_location( const Vec3D &camera )
{
    _camera = camera;
}


void PerspectiveCamera::set_camera_up( const Vec3D &up )
{
    _up = up;
}


void PerspectiveCamera::set_zplanes( double near, double far )
{
    _near = near;
    _far = far;
}



