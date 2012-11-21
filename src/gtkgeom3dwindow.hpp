/*! \file gtkgeom3dwindow.hpp
 *  \brief %Geometry view window for 3d
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


#ifndef GTKGEOM3DWINDOW_HPP
#define GTKGEOM3DWINDOW_HPP 1


#include <gtk/gtk.h>
#include "gtkwindow.hpp"
#include "gtkplotter.hpp"
#include "transformation.hpp"


/*! \brief Interactive geometry plotter window.
 */
class GTKGeom3DWindow : public GTKWindow {

    GTKPlotter             &_plotter;
    const Geometry         &_geom;
    const ParticleDataBase *_pdb;

    GtkWidget              *_window;
    GtkWidget              *_darea;
    GtkWidget              *_menubar;
    GtkWidget              *_menu_file;
    GtkWidget              *_toolbar;
    GtkWidget              *_statusbar;

    GtkToolItem            *_radioitem;

    size_t                  _width;
    size_t                  _height;

    int                     _tool;
    double                  _oldx;
    double                  _oldy;
    double                  _endx;
    double                  _endy;

    double                  _near;
    double                  _far;
    Vec3D                   _camera;
    Vec3D                   _target;
    Vec3D                   _up;
    double                  _zoom;
    Transformation          _modeltrans;

    class Renderer         *_renderer;
    Vec3D                   _center;
    double                  _scale;

    int32_t                 _clevel[6];
    std::vector<Vec3D>      _csurface[6];  // Sequences of 3: x0,x1,x2
    std::vector<Vec3D>      _gsurface;     // Sequences of 4: norm,x0,x1,x2

    int32_t                 _pdiv;
    bool                    _bbox;

    void init_window( void );
    void init_renderer( void );
    void init_model( void );
    void init_camera_and_rotation( void );

    void cplane_add_vertex( int32_t p, const int32_t i[3], const int32_t vb[3], double dx, double dy );
    double cplane_dist( const int32_t i[3], const int32_t vb[3], int32_t dx, int32_t dy, int32_t dir );
    int32_t case2d( const int32_t i[3], const int32_t vb[3] );
    void build_cut_plane( int32_t p, const int32_t vb[3], int32_t level );
    void build_cut_planes( void );
    void build_geometry_surface( void );
    void clear_surface_data( void );

    void draw_cut_planes( void );
    void draw_model( void );
    void draw_bbox( void );
    void draw_beam( void );

    void setup_lights( void );

    void move( int action, double x, double y );
    void zoom_out( double x, double y );
    void zoom_in( double x, double y );
    void zoom_window( int action, double x, double y );
    void zoom_fit( void );
    void geom2d_launch( void );
    void hardcopy( void );


    void track( int action, double x, double y );

    void darea_motion( GdkEventMotion *event );
    void darea_enter( GdkEventCrossing *event );
    void darea_leave( GdkEventCrossing *event );
    void darea_button( GdkEventButton *event );

    void expose( void );
    void configure( void );
    void delete_window( void );

    void menuitem_tool_change( GtkToolButton *button );
    void menuitem_preferences( GtkMenuItem *menuitem );

    static gboolean window_delete_signal( GtkWidget *widget, 
					  GdkEventExpose *event, 
					  gpointer object );

    static void menuitem_quit_signal( GtkMenuItem *menuitem,
				      gpointer object );
    static void menuitem_tool_change_signal( GtkToolButton *button,
					     gpointer object );
    static void menuitem_hardcopy_signal( GtkToolButton *button,
					  gpointer object );
    static void menuitem_zoom_fit_signal( GtkToolButton *button,
					  gpointer object );
    static void menuitem_geom2d_signal( GtkToolButton *button,
					gpointer object );
    static void menuitem_preferences_signal( GtkMenuItem *menuitem,
					     gpointer object );

    static gboolean darea_configure_signal( GtkWidget *widget, 
					    GdkEventConfigure *event, 
					    gpointer object );
    static gboolean darea_expose_signal( GtkWidget *widget, 
					 GdkEventExpose *event, 
					 gpointer object );
    static gboolean darea_button_signal( GtkWidget *widget, 
					 GdkEventButton *event,
					 gpointer object );
    static gboolean darea_motion_signal( GtkWidget *widget, 
					 GdkEventMotion *event,
					 gpointer object );
    static gboolean darea_enter_signal( GtkWidget *widget, 
					GdkEventCrossing *event,
					gpointer object );
    static gboolean darea_leave_signal( GtkWidget *widget, 
					GdkEventCrossing *event,
					gpointer object );

public:

    /*! \brief Constructor.
     */
    GTKGeom3DWindow( GTKPlotter &plotter,
		     const Geometry &geom,
		     const ParticleDataBase *pdb );

    /*! \brief Destructor.
     */
    virtual ~GTKGeom3DWindow();
};


#endif
