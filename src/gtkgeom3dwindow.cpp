/*! \file gtkgeom3dwindow.cpp
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


#include <cairo.h>
#include "config.h"

#ifdef OPENGL
#include "glrenderer.hpp"
#endif

#include "softwarerenderer.hpp"
#include "gtkgeom3dwindow.hpp"
#include "gtkhardcopy.hpp"
#include "icons.hpp"
#include "ibsimu.hpp"


#define TOOL_UNKNOWN  -1
#define TOOL_MOVE      0
#define TOOL_ZOOM_IN   1
#define TOOL_ZOOM_OUT  2
#define TOOL_TRACK     3


GTKGeom3DWindow::GTKGeom3DWindow( GTKPlotter &plotter,
				  const Geometry &geom,
				  const ParticleDataBase *pdb )
    : _plotter(plotter), _geom3dplot(geom,pdb), _geom(geom), _width(640), _height(480)
{
    if( !geom.surface_built() ) {
	throw( Error( ERROR_LOCATION, "geometry surface not built" ) );	
    }

    init_window();
    init_renderer();

    g_signal_connect( G_OBJECT(_darea), "configure_event",
		      G_CALLBACK(darea_configure_signal), 
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "expose_event",
		      G_CALLBACK(darea_expose_signal), 
		      (gpointer)this );

    gtk_widget_show_all( _window );
}


void GTKGeom3DWindow::init_renderer( void )
{
#ifdef OPENGL
    if( _plotter.opengl() ) {
	try {
	    // Try initializing OpenGL renderer
	    _renderer = new GLRenderer( _darea );
	} catch( GLRenderer::ErrorGLInit e ) {
	    // Fallback to software renderer
	    _renderer = new SoftwareRenderer( _darea );
	}
    } else {
	// No GdkGLExt initialized
	_renderer = new SoftwareRenderer( _darea );
    }
#else
    // Only software renderer available
    _renderer = new SoftwareRenderer( _darea );
#endif
}


void GTKGeom3DWindow::init_window( void )
{
    // Window
    _window = gtk_window_new( GTK_WINDOW_TOPLEVEL );
    gtk_window_set_title( GTK_WINDOW(_window), "Simulation 3d geometry" );
    g_signal_connect( G_OBJECT(_window), "delete_event",
		      G_CALLBACK(window_delete_signal), 
		      (gpointer)this );
    GtkWidget *vbox;
    vbox = gtk_vbox_new( FALSE, 0 );

    // Menu bar
    _menubar = gtk_menu_bar_new();
    _menu_file = gtk_menu_new();
    GtkWidget *item_hardcopy = gtk_menu_item_new_with_mnemonic( "_Hardcopy" );
    gtk_menu_shell_append( GTK_MENU_SHELL(_menu_file), item_hardcopy );
    GtkWidget *item_quit = gtk_menu_item_new_with_mnemonic( "_Quit" );
    gtk_menu_shell_append( GTK_MENU_SHELL(_menu_file), item_quit );
    GtkWidget *item_file = gtk_menu_item_new_with_mnemonic( "_File" );
    gtk_menu_item_set_submenu( GTK_MENU_ITEM(item_file), _menu_file );
    gtk_menu_bar_append( GTK_MENU_BAR(_menubar), item_file );
    g_signal_connect( G_OBJECT(item_quit), "activate",
		      G_CALLBACK(menuitem_quit_signal),
		      (gpointer)this );
    g_signal_connect( G_OBJECT(item_hardcopy), "activate",
		      G_CALLBACK(menuitem_hardcopy_signal),
		      (gpointer)this );
    gtk_box_pack_start( GTK_BOX(vbox), _menubar, FALSE, TRUE, 0 );

    // Add Edit/Configure menu
    GtkWidget *menu_edit, *item_edit, *item_preferences;
    menu_edit = gtk_menu_new();
    item_preferences = gtk_menu_item_new_with_mnemonic( "_Preferences" );
    gtk_menu_shell_append( GTK_MENU_SHELL(menu_edit), item_preferences );
    item_edit = gtk_menu_item_new_with_mnemonic( "_Edit" );
    gtk_menu_item_set_submenu( GTK_MENU_ITEM(item_edit), menu_edit );
    gtk_menu_shell_insert( GTK_MENU_SHELL(_menubar), item_edit, 1 );
    g_signal_connect( G_OBJECT(item_preferences), "activate",
                      G_CALLBACK(menuitem_preferences_signal),
                      (gpointer)this );

    // Tool bar
    _toolbar = gtk_toolbar_new();
    gtk_toolbar_set_style( GTK_TOOLBAR(_toolbar), GTK_TOOLBAR_ICONS );
    gtk_toolbar_set_orientation( GTK_TOOLBAR(_toolbar), GTK_ORIENTATION_HORIZONTAL );
    

    // Creating "Hardcopy" button
    GdkPixbuf *pixbuf = gdk_pixbuf_new_from_inline( -1, icon_hardcopy_inline, FALSE, NULL );
    GtkWidget *icon = gtk_image_new_from_pixbuf( pixbuf );
    GtkToolItem *toolitem = gtk_tool_button_new( icon, "Hardcopy" );
#if GTK_CHECK_VERSION(2,12,0)
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Hardcopy" );
#endif
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "clicked",
		      G_CALLBACK(menuitem_hardcopy_signal),
		      (gpointer)this );
    
    // Creating separator
    toolitem = gtk_separator_tool_item_new();
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    
    // Creating "Zoom in" button
    pixbuf = gdk_pixbuf_new_from_inline( -1, icon_zoom_in_inline, FALSE, NULL );
    icon = gtk_image_new_from_pixbuf( pixbuf );
    _radioitem = toolitem = gtk_radio_tool_button_new( NULL );
    gtk_tool_button_set_label( GTK_TOOL_BUTTON(toolitem), "Zoom in" );
#if GTK_CHECK_VERSION(2,12,0)
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Zoom in" );
#endif
    gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
    gtk_toggle_tool_button_set_active( GTK_TOGGLE_TOOL_BUTTON(toolitem), FALSE );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "toggled",
		      G_CALLBACK(menuitem_tool_change_signal),
		      (gpointer)this );

    // Creating "Zoom out" button
    pixbuf = gdk_pixbuf_new_from_inline( -1, icon_zoom_out_inline, FALSE, NULL );
    icon = gtk_image_new_from_pixbuf( pixbuf );
    toolitem = gtk_radio_tool_button_new_from_widget( GTK_RADIO_TOOL_BUTTON(_radioitem) );
    gtk_tool_button_set_label( GTK_TOOL_BUTTON(toolitem), "Zoom out" );
#if GTK_CHECK_VERSION(2,12,0)
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Zoom out" );
#endif
    gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
    gtk_toggle_tool_button_set_active( GTK_TOGGLE_TOOL_BUTTON(toolitem), FALSE );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "toggled",
		      G_CALLBACK(menuitem_tool_change_signal),
		      (gpointer)this );

    // Creating "Zoom fit" button
    pixbuf = gdk_pixbuf_new_from_inline( -1, icon_zoom_fit_inline, FALSE, NULL );
    icon = gtk_image_new_from_pixbuf( pixbuf );
    toolitem = gtk_tool_button_new( icon, "Zoom fit" );
#if GTK_CHECK_VERSION(2,12,0)
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Zoom fit" );
#endif
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "clicked",
		      G_CALLBACK(menuitem_zoom_fit_signal),
		      (gpointer)this );
    
    // Creating "Move" button
    pixbuf = gdk_pixbuf_new_from_inline( -1, icon_move_inline, FALSE, NULL );
    icon = gtk_image_new_from_pixbuf( pixbuf );
    toolitem = gtk_radio_tool_button_new_from_widget( GTK_RADIO_TOOL_BUTTON(_radioitem) );
    gtk_tool_button_set_label( GTK_TOOL_BUTTON(toolitem), "Move" );
#if GTK_CHECK_VERSION(2,12,0)
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Move" );
#endif
    gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
    gtk_toggle_tool_button_set_active( GTK_TOGGLE_TOOL_BUTTON(toolitem), TRUE );
    _tool = TOOL_MOVE;
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "toggled",
		      G_CALLBACK(menuitem_tool_change_signal),
		      (gpointer)this );

    // Creating separator
    toolitem = gtk_separator_tool_item_new();
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    
    // Creating "Track" button
    pixbuf = gdk_pixbuf_new_from_inline( -1, icon_track_inline, FALSE, NULL );
    icon = gtk_image_new_from_pixbuf( pixbuf );
    toolitem = gtk_radio_tool_button_new_from_widget( GTK_RADIO_TOOL_BUTTON(_radioitem) );
    gtk_tool_button_set_label( GTK_TOOL_BUTTON(toolitem), "Track" );
#if GTK_CHECK_VERSION(2,12,0)
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Track plot" );
#endif
    gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "toggled",
		      G_CALLBACK(menuitem_tool_change_signal),
		      (gpointer)this );
    
    // Creating separator
    toolitem = gtk_separator_tool_item_new();
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );

    // Creating "Geom 2D" button
    pixbuf = gdk_pixbuf_new_from_inline( -1, icon_geom2d_inline, FALSE, NULL );
    icon = gtk_image_new_from_pixbuf( pixbuf );
    toolitem = gtk_tool_button_new( icon, "2D geometry view" );
#if GTK_CHECK_VERSION(2,12,0)
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "2D geometry view" );
#endif
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "clicked",
		      G_CALLBACK(menuitem_geom2d_signal),
		      (gpointer)this );

    // Add toolbar to vbox
    gtk_box_pack_start( GTK_BOX(vbox), _toolbar, FALSE, TRUE, 0 );

    // Statusbar
    _statusbar = gtk_statusbar_new();
    gtk_statusbar_push( GTK_STATUSBAR(_statusbar), 0, "Done" );
    gtk_box_pack_end( GTK_BOX(vbox), _statusbar, FALSE, TRUE, 0 );


    // Drawing area
    _darea = gtk_drawing_area_new();
    gtk_widget_set_size_request( _darea, _width, _height );
    gtk_widget_add_events( _darea, GDK_EXPOSURE_MASK |
			   GDK_LEAVE_NOTIFY_MASK |
			   GDK_ENTER_NOTIFY_MASK |
			   GDK_POINTER_MOTION_HINT_MASK |
			   GDK_POINTER_MOTION_MASK |
			   GDK_SCROLL_MASK |
			   GDK_BUTTON_PRESS_MASK |
			   GDK_BUTTON_RELEASE_MASK |
			   GDK_BUTTON_MOTION_MASK );

    g_signal_connect( G_OBJECT(_darea), "button_press_event",
		      G_CALLBACK(darea_button_signal),
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "button_release_event",
		      G_CALLBACK(darea_button_signal),
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "motion_notify_event",
		      G_CALLBACK(darea_motion_signal),
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "enter_notify_event",
		      G_CALLBACK(darea_enter_signal),
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "leave_notify_event",
		      G_CALLBACK(darea_leave_signal),
		      (gpointer)this );
    gtk_box_pack_end( GTK_BOX(vbox), _darea, TRUE, TRUE, 0 );

    // Add vbox to window
    gtk_container_add( GTK_CONTAINER(_window), vbox );

    gtk_window_present( GTK_WINDOW(_window) );
}


GTKGeom3DWindow::~GTKGeom3DWindow()
{
}


void GTKGeom3DWindow::delete_window( void )
{
    _plotter.delete_window( this );
}


void GTKGeom3DWindow::configure( void )
{
    GtkAllocation alloc;
    gtk_widget_get_allocation( _darea, &alloc );
    _width = alloc.width;
    _height = alloc.height;
    _geom3dplot.set_size( _width, _height );
}


void GTKGeom3DWindow::expose( void )
{
    _geom3dplot.draw( _renderer );
}


// action 0: button press (move start), 1: pointer move, 
// 2: release button (move end)
void GTKGeom3DWindow::move( int action, double x, double y )
{
    if( action == 0 ) {
	_oldx = x;
	_oldy = y;
    } else {
	double dx = x-_oldx;
	double dy = y-_oldy;
	Transformation modeltrans = _geom3dplot.get_model_transformation();
	modeltrans.rotate_y( 3.0*dx/_width );
	modeltrans.rotate_x( 3.0*dy/_height );
	_geom3dplot.set_model_transformation( modeltrans );
	_oldx = x;
	_oldy = y;
	gtk_widget_queue_draw_area( _darea, 0, 0, _width, _height );
    }
}


void GTKGeom3DWindow::zoom_out( double x, double y )
{
    double near, far, zoom;
    _geom3dplot.get_projection_frustum( near, far, zoom );
    Vec3D target, camera, up;
    _geom3dplot.get_view_look_at( camera, target, up );

    double fac = 1.414;
    double left = -zoom*_width/_height;
    double right= zoom*_width/_height;
    double top = zoom;
    double bottom = -zoom;
    double u = (left + (right-left)*(x/_width));
    double v = (bottom + (top-bottom)*(1.0-y/_height));
    Vec3D viewdir = target-camera;
    viewdir.normalize();
    viewdir *= near;
    Vec3D rightdir = cross(viewdir,up);
    Vec3D real_up = cross(rightdir,viewdir);
    rightdir.normalize();
    real_up.normalize();
    target = camera + near*viewdir + u*rightdir + v*real_up;
    zoom *= fac;

    _geom3dplot.set_projection_zoom( zoom );
    _geom3dplot.set_view_look_at( camera, target, up );
    gtk_widget_queue_draw_area( _darea, 0, 0, _width, _height );
}


void GTKGeom3DWindow::zoom_in( double x, double y )
{
    double near, far, zoom;
    _geom3dplot.get_projection_frustum( near, far, zoom );
    Vec3D target, camera, up;
    _geom3dplot.get_view_look_at( camera, target, up );

    double fac = 0.707;
    double left = -zoom*_width/_height;
    double right= zoom*_width/_height;
    double top = zoom;
    double bottom = -zoom;
    double u = (left + (right-left)*(x/_width));
    double v = (bottom + (top-bottom)*(1.0-y/_height));
    Vec3D viewdir = target-camera;
    viewdir.normalize();
    viewdir *= near;
    Vec3D rightdir = cross(viewdir,up);
    Vec3D real_up = cross(rightdir,viewdir);
    rightdir.normalize();
    real_up.normalize();
    target = camera + near*viewdir + u*rightdir + v*real_up;
    zoom *= fac;

    _geom3dplot.set_projection_zoom( zoom );
    _geom3dplot.set_view_look_at( camera, target, up );
    gtk_widget_queue_draw_area( _darea, 0, 0, _width, _height );
}


// action: 0 for press, 1 for modify, 2 for release
void GTKGeom3DWindow::zoom_window( int action, double x, double y )
{
    if( action == 0 ) {
	//std::cout << "Zoom window start\n";
	_oldx = _endx = (int)floor(x+0.5);
	_oldy = _endy = (int)floor(y+0.5);
    } else if( action == 1 ) {
	//std::cout << "Zoom window modify\n";
    } else {
	//std::cout << "Zoom window end\n";
    }
}


void GTKGeom3DWindow::zoom_fit( void )
{
    _geom3dplot.reset_camera_and_rotation();
    gtk_widget_queue_draw_area( _darea, 0, 0, _width, _height );
}


// action: 0 for leave, 1 for motion, 2 for enter, 3 for create, 4 for delete
void GTKGeom3DWindow::track( int action, double x, double y )
{

}


void GTKGeom3DWindow::darea_motion( GdkEventMotion *event )
{
    if( (_tool == TOOL_MOVE && (event->state & GDK_BUTTON1_MASK)) || 
	event->state & GDK_BUTTON2_MASK ) {
	move( 1, event->x, event->y );
    } else if( _tool == TOOL_TRACK ) {
	track( 1, event->x, event->y );
    } else if( (_tool == TOOL_ZOOM_IN && (event->state & GDK_BUTTON1_MASK)) ||
	       (event->state & GDK_BUTTON2_MASK) ) {
	zoom_window( 1, event->x, event->y );
    }
}


void GTKGeom3DWindow::darea_enter( GdkEventCrossing *event )
{
    if( _tool == TOOL_TRACK )
	track( 2, 0, 0 );
}


void GTKGeom3DWindow::darea_leave( GdkEventCrossing *event )
{
    if( _tool == TOOL_TRACK )
	track( 0, event->x, event->y );
}


void GTKGeom3DWindow::darea_button( GdkEventButton *event )
{
    if( (_tool == TOOL_MOVE && event->type == GDK_BUTTON_PRESS && event->button == 1) || 
	(event->type == GDK_BUTTON_PRESS && event->button == 2)	) {
	move( 0, event->x, event->y );
    } else if( (_tool == TOOL_MOVE && event->type == GDK_BUTTON_RELEASE && event->button == 1) || 
	       (event->type == GDK_BUTTON_RELEASE && event->button == 2) ) {
	move( 2, event->x, event->y );
    } else if( _tool == TOOL_ZOOM_OUT && event->type == GDK_BUTTON_RELEASE && event->button == 1 ) {
	zoom_out( event->x, event->y );
    } else if( _tool == TOOL_ZOOM_IN && event->type == GDK_BUTTON_PRESS && event->button == 1 ) {
	zoom_window( 0, event->x, event->y );
    } else if( _tool == TOOL_ZOOM_IN && event->type == GDK_BUTTON_RELEASE && event->button == 1 ) {
	if( fabs( _oldx - event->x ) < 2.0 ||
	    fabs( _oldy - event->y ) < 2.0 )
	    zoom_in( event->x, event->y );
	else
	    zoom_window( 2, event->x, event->y );
    }
}


void GTKGeom3DWindow::menuitem_tool_change( GtkToolButton *button )
{
    int tool;
    const char *label = gtk_tool_button_get_label( button );
    if( !strcmp( label, "Zoom in" ) ) {
	tool = TOOL_ZOOM_IN;
    } else if( !strcmp( label, "Zoom out" ) ) {
	tool = TOOL_ZOOM_OUT;
    } else if( !strcmp( label, "Move" ) ) {
	tool = TOOL_MOVE;
    } else if( !strcmp( label, "Track" ) ) {
	tool = TOOL_TRACK;
    } else {
	tool = TOOL_UNKNOWN;
	return;
    }

    if( !gtk_toggle_tool_button_get_active( GTK_TOGGLE_TOOL_BUTTON(button) ) ) {
	// Disable tool
	if( tool == TOOL_TRACK )
	    track( 4, 0, 0 );
	_tool = TOOL_UNKNOWN;
    } else {
	// Enable tool
	_tool = tool;
	if( tool == TOOL_TRACK )
	    track( 3, 0, 0 );
    }
}


void GTKGeom3DWindow::menuitem_preferences( GtkMenuItem *menuitem )
{
    GtkWidget *dialog = gtk_dialog_new_with_buttons( "Plot preferences",
						     GTK_WINDOW(_window),
						     (GtkDialogFlags)(GTK_DIALOG_MODAL | 
								      GTK_DIALOG_DESTROY_WITH_PARENT), 
						     GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
						     GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
						     NULL );
    GtkWidget *vbox = GTK_DIALOG(dialog)->vbox;
    GtkWidget *notebook = gtk_notebook_new();

    // ****************************************************************************

    GtkWidget *vbox2 = gtk_vbox_new( FALSE, 0 );

    GtkWidget *label = gtk_label_new( "Trajectory division" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    GtkWidget *spinbutton_pdiv = gtk_spin_button_new_with_range( 0, 1000000, 1.0 );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(spinbutton_pdiv), _geom3dplot.get_particle_div() );
    GtkWidget *hbox = gtk_hbox_new( TRUE, 30 );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), spinbutton_pdiv, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox2), hbox, FALSE, TRUE, 0 );

    label = gtk_label_new( "Trajectory offset" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    GtkWidget *spinbutton_poffset = gtk_spin_button_new_with_range( 0, 1000000, 1.0 );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(spinbutton_poffset), _geom3dplot.get_particle_offset() );
    hbox = gtk_hbox_new( TRUE, 30 );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), spinbutton_poffset, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox2), hbox, FALSE, TRUE, 0 );

    label = gtk_label_new( "BBox" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    GtkWidget *button_bbox = gtk_check_button_new_with_label( "on/off" );
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(button_bbox), _geom3dplot.get_bbox() );
    hbox = gtk_hbox_new( TRUE, 30 );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), button_bbox, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox2), hbox, FALSE, TRUE, 0 );

    // Notebook page
    label = gtk_label_new( "Misc" );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );

    // ****************************************************************************

    vbox2 = gtk_vbox_new( FALSE, 0 );

    // Cut levels
    GtkWidget *spinbutton[6];
    for( int a = 0; a < 6; a++ ) {
	hbox = gtk_hbox_new( TRUE, 30 );
	if( a == 0 ) {
	    label = gtk_label_new( "Cut level xmin" );
	    spinbutton[a] = gtk_spin_button_new_with_range( 0, _geom.size(0)-1, 1.0 );
	} else if( a == 1 ) {
	    label = gtk_label_new( "Cut level xmax" );
	    spinbutton[a] = gtk_spin_button_new_with_range( 0, _geom.size(0)-1, 1.0 );
	} else if( a == 2 ) {
	    label = gtk_label_new( "Cut level ymin" );
	    spinbutton[a] = gtk_spin_button_new_with_range( 0, _geom.size(1)-1, 1.0 );
	} else if( a == 3 ) {
	    label = gtk_label_new( "Cut level ymax" );
	    spinbutton[a] = gtk_spin_button_new_with_range( 0, _geom.size(1)-1, 1.0 );
	} else if( a == 4 ) {
	    label = gtk_label_new( "Cut level zmin" );
	    spinbutton[a] = gtk_spin_button_new_with_range( 0, _geom.size(2)-1, 1.0 );
	} else if( a == 5 ) {
	    label = gtk_label_new( "Cut level zmax" );
	    spinbutton[a] = gtk_spin_button_new_with_range( 0, _geom.size(2)-1, 1.0 );
	}
	gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON(spinbutton[a]), _geom3dplot.get_clevel( a ) );
	gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
	gtk_box_pack_start( GTK_BOX(hbox), spinbutton[a], FALSE, TRUE, 0 );
	gtk_box_pack_start( GTK_BOX(vbox2), hbox, FALSE, TRUE, 0 );
    }

    // Notebook page
    label = gtk_label_new( "Cut levels" );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );

    // ****************************************************************************

    // Pack notebook
    gtk_box_pack_start( GTK_BOX(vbox), notebook, FALSE, TRUE, 0 );

    gtk_widget_show_all( dialog );
    if( gtk_dialog_run( GTK_DIALOG(dialog) ) == GTK_RESPONSE_ACCEPT ) {

	// Misc
	uint32_t particle_div = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(spinbutton_pdiv) );
	uint32_t particle_offset = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(spinbutton_poffset) );
	_geom3dplot.set_particle_div( particle_div, particle_offset );
	_geom3dplot.set_bbox( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(button_bbox) ) );

	// Cut levels
	for( int a = 0; a < 6; a++ )
	    _geom3dplot.set_clevel( a, gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(spinbutton[a]) ) );

	// Rebuild model and refresh output
	_geom3dplot.rebuild_model();
	gtk_widget_queue_draw_area( _darea, 0, 0, _width, _height );
    }

    gtk_widget_destroy( dialog );
}


void GTKGeom3DWindow::geom2d_launch( void )
{
    _plotter.new_geometry_plot_window();
}


void GTKGeom3DWindow::hardcopy( void )
{
    GTKHardcopy hardcopy( _window, &_geom3dplot, _width, _height );
    hardcopy.run();
}


/* **********************************************
 * STATIC SIGNAL FUNCTIONS 
 */

gboolean GTKGeom3DWindow::darea_configure_signal( GtkWidget *widget, 
						  GdkEventConfigure *event, 
						  gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->configure();
    return( FALSE );
}


gboolean GTKGeom3DWindow::darea_expose_signal( GtkWidget *widget, 
					       GdkEventExpose *event, 
					       gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->expose();
    return( FALSE );
}


gboolean GTKGeom3DWindow::darea_motion_signal( GtkWidget *widget, 
					       GdkEventMotion *event,
					       gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->darea_motion( event );
    return( FALSE );  
}


gboolean GTKGeom3DWindow::darea_enter_signal( GtkWidget *widget, 
					      GdkEventCrossing *event,
					      gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->darea_enter( event );
    return( FALSE );  
}


gboolean GTKGeom3DWindow::darea_leave_signal( GtkWidget *widget, 
					      GdkEventCrossing *event,
					      gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->darea_leave( event );
    return( FALSE );
}


gboolean GTKGeom3DWindow::darea_button_signal( GtkWidget *widget, 
					       GdkEventButton *event,
					       gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->darea_button( event );
    return( FALSE );
}


gboolean GTKGeom3DWindow::window_delete_signal( GtkWidget *widget, 
						GdkEventExpose *event, 
						gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->delete_window();
    return( FALSE );
}


void GTKGeom3DWindow::menuitem_preferences_signal( GtkMenuItem *menuitem,
						   gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->menuitem_preferences( menuitem );
}


void GTKGeom3DWindow::menuitem_quit_signal( GtkMenuItem *menuitem,
					    gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->delete_window();
}


void GTKGeom3DWindow::menuitem_tool_change_signal( GtkToolButton *button,
						   gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->menuitem_tool_change( button );
}


void GTKGeom3DWindow::menuitem_geom2d_signal( GtkToolButton *button,
					    gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->geom2d_launch();
}


void GTKGeom3DWindow::menuitem_hardcopy_signal( GtkToolButton *button,
						gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->hardcopy();
}


void GTKGeom3DWindow::menuitem_zoom_fit_signal( GtkToolButton *button,
						gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->zoom_fit();
}

