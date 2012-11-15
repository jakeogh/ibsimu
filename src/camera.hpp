/*! \file camera.hpp
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


#ifndef CAMERA_HPP
#define CAMERA_HPP 1


#include "vec3d.hpp"
#include "transformation.hpp"


/*! \brief Camera.
 *
 */
class Camera {
    
public:
    
    /*! \brief Destructor.
     */
    virtual ~Camera() {}
    
    /*! \brief Initialize camera.
     *
     *  Sets up viewport, projection and modelview matrices for
     *  camera view of the geometry in OpenGL.
     */
    virtual void gl_initalize_camera( void ) = 0;
    
    /*! \brief Set up viewport size.
     */
    virtual void set_size( double width, double height ) = 0;
    
    /*! \brief Change view relative to current.
     *
     *  Sets \a (x,y) as the new center of view and changes zoom
     *  setting with factor \a fac.
     */
    virtual void set_view_relative( double x, double y, double fac ) = 0;

    /*! \brief Change field of view.
     *
     *  The field of view is set to \a +/-zoom at distance \a near
     *  from the camera in the longer direction (horizontal/vertical).
     */
    virtual void set_field_of_view( double zoom ) = 0;
    
    /*! \brief Change view target.
     */
    virtual void set_target_location( const Vec3D &target ) = 0;
    
    /*! \brief Change camera location.
     */
    virtual void set_camera_location( const Vec3D &camera ) = 0;
    
    /*! \brief Change camera up direction.
     */
    virtual void set_camera_up( const Vec3D &up ) = 0;
    
    /*! \brief Change near and far planes.
     */
    virtual void set_zplanes( double near, double far ) = 0;
};


/*! \brief Constructor.
 */
class PerspectiveCamera : public Camera {
    
    double          _width;
    double          _height;
    double          _zoom;
    double          _near;
    double          _far;
    double          _left;
    double          _right;
    double          _bottom;
    double          _top;
    Vec3D           _target;
    Vec3D           _up;
    Vec3D           _camera;
    
public:
    
    /*! \brief Constructor.
     */
    PerspectiveCamera();
    
    /*! \brief Destructor.
     */
    virtual ~PerspectiveCamera();
    
    /*! \brief Initialize camera.
     *
     *  Sets up viewport, projection and modelview matrices for
     *  camera view of the geometry in OpenGL.
     */
    virtual void gl_initalize_camera( void );
    
    /*! \brief Set up viewport size.
     */
    virtual void set_size( double width, double height );
    
    /*! \brief Change view relative to current.
     *
     *  Sets \a (x,y) as the new center of view and changes zoom
     *  setting by multiplying with factor \a fac.
     */
    virtual void set_view_relative( double x, double y, double fac );

    /*! \brief Change field of view.
     *
     *  The field of view is set to \a +/-zoom at distance \a near
     *  from the camera in the longer direction (horizontal/vertical).
     */
    virtual void set_field_of_view( double zoom );
    
    /*! \brief Change view target.
     */
    virtual void set_target_location( const Vec3D &target );
    
    /*! \brief Change camera location.
     */
    virtual void set_camera_location( const Vec3D &camera );
    
    /*! \brief Change camera up direction.
     */
    virtual void set_camera_up( const Vec3D &up );
    
    /*! \brief Change near and far planes.
     */
    virtual void set_zplanes( double near, double far );
};


#endif
