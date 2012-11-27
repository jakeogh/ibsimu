/*! \file geom3dplot.hpp
 *  \brief %Geometry 3d plotter
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


#ifndef GEOM3DPLOT_HPP
#define GEOM3DPLOT_HPP 1


#include <stdint.h>
#include <vector>
#include "geometry.hpp"
#include "particledatabase.hpp"
#include "transformation.hpp"
#include "vec3d.hpp"


/*! \brief %Geometry 3D plotter.
 */
class Geom3DPlot {

    const Geometry         &_geom;
    const ParticleDataBase *_pdb;

    size_t                  _width;
    size_t                  _height;

    double                  _near;
    double                  _far;
    Vec3D                   _camera;
    Vec3D                   _target;
    Vec3D                   _up;
    double                  _zoom;

    Transformation          _modeltrans;

    Vec3D                   _center;
    double                  _scale;

    std::vector<Vec3D>      _csurface[6];  // Sequences of 3: x0,x1,x2
    std::vector<Vec3D>      _gsurface;     // Sequences of 4: norm,x0,x1,x2

    uint32_t                _particle_div;
    uint32_t                _particle_offset;
    bool                    _bbox;
    uint32_t                _clevel[6];

    void cplane_add_vertex( int32_t p, const int32_t i[3], const int32_t vb[3], double dx, double dy );
    double cplane_dist( const int32_t i[3], const int32_t vb[3], int32_t dx, int32_t dy, int32_t dir );
    int32_t case2d( const int32_t i[3], const int32_t vb[3] );
    void build_cut_plane( int32_t p, const int32_t vb[3], int32_t level );
    void build_cut_planes( void );
    void build_geometry_surface( void );
    void clear_surface_data( void );

    void draw_cut_planes( class Renderer *r );
    void draw_model( class Renderer *r );
    void draw_bbox( class Renderer *r );
    void draw_beam( class Renderer *r );

public:

    /*! \brief Constructor.
     */
    Geom3DPlot( const Geometry &geom,
		const ParticleDataBase *pdb );

    /*! \brief Destructor.
     */
    ~Geom3DPlot();

    /*! \brief Set canvas size.
     */
    void set_size( uint32_t width, uint32_t height );

    /*! \brief Set modelling transformation.
     */
    void set_model_transformation( const Transformation &modeltrans );

    /*! \brief Get modelling transformation.
     */
    Transformation get_model_transformation( void ) const;

    /*! \brief Set camera projection.
     */
    void set_projection_frustum( double near, double far, double zoom );
    
    /*! \brief Get camera projection.
     */
    void get_projection_frustum( double &near, double &far, double &zoom ) const ;

    /*! \brief Set camera zoom
     */
    void set_projection_zoom( double zoom );

    /*! \brief Get camera zoom
     */
    double get_projection_zoom( void ) const;

    /*! \brief Set camera positioning.
     */
    void set_view_look_at( const Vec3D &camera, 
			   const Vec3D &target,
			   const Vec3D &up );

    /*! \brief Get camera positioning.
     */
    void get_view_look_at( Vec3D &camera, 
			   Vec3D &target,
			   Vec3D &up ) const;

    /*! \brief Rebuild model.
     *
     *  Needed after changing settings affecting the displayed model.
     */
    void rebuild_model( void );

    /*! \brief Reset camera settings and modelling transformation to default.
     */
    void reset_camera_and_rotation( void );

    /*! \brief Set particle divisor and offset.
     *
     *  Set \a particle_div to zero for no plotting, one for plotting
     *  every particle, two for plotting every second particle, three
     *  for plotting every third particle, etc. Defaults to
     *  11. Plotter skips the first \a particle_offset particles.
     */
    void set_particle_div( uint32_t particle_div, uint32_t particle_offset = 0 );
    
    /*! \brief Get particle division.
     */
    uint32_t get_particle_div( void ) const;

    /*! \brief Get particle offset.
     */
    uint32_t get_particle_offset( void ) const;

    /*! \brief Set if bounding box is plotted.
     */
    void set_bbox( bool bbox );

    /*! \brief Get if bounding box is plotted.
     */
    bool get_bbox( void ) const;

    /*! \brief Set cut level for direction \a a.
     *
     *  Direction 0 is xmin, direction 1 is xmax, direction 2 is ymin, direction 3 is
     *  ymax, direction 4 is zmin and direction 5 is zmax.
     */
    void set_clevel( uint32_t direction, uint32_t level );

    /*! \brief Get cut level for direction \a a.
     *
     *  Direction 0 is xmin, direction 1 is xmax, direction 2 is ymin, direction 3 is
     *  ymax, direction 4 is zmin and direction 5 is zmax.
     */
    uint32_t get_clevel( uint32_t direction ) const;

    /*! \brief Draw 3D geometry using renderer \a r.
     */
    void draw( class Renderer *r );
};


#endif
