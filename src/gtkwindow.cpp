/*! \file gtkwindow.cpp
 *  \brief Source code for gtkwindow.cpp
 */

/* Copyright (c) 2005-2009 Taneli Kalvas. All rights reserved.
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
 * tvkalvas@cc.jyu.fi.
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
#include "gtkwindow.hpp"
#include "gtkhardcopy.hpp"
#include "gtkpreferences.hpp"
#include "icons.hpp"



#define TOOL_UNKNOWN  -1
#define TOOL_MOVE      0
#define TOOL_ZOOM_IN   1
#define TOOL_ZOOM_OUT  2
#define TOOL_TRACK     3



GTKWindow::GTKWindow( class GTKPlotter *plotter )
    : _width(640), _height(480), _cairo(NULL), _surface(NULL), _plotter(plotter)
{
    //std::cout << "GTKWindow constructor\n";

    // Window
    _window = gtk_window_new( GTK_WINDOW_TOPLEVEL );
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
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Hardcopy" );
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
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Zoom in" );
    gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
    gtk_toggle_tool_button_set_active( GTK_TOGGLE_TOOL_BUTTON(toolitem), TRUE );
    _tool = TOOL_ZOOM_IN;
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "toggled",
		      G_CALLBACK(menuitem_tool_change_signal),
		      (gpointer)this );

    // Creating "Zoom out" button
    pixbuf = gdk_pixbuf_new_from_inline( -1, icon_zoom_out_inline, FALSE, NULL );
    icon = gtk_image_new_from_pixbuf( pixbuf );
    toolitem = gtk_radio_tool_button_new_from_widget( GTK_RADIO_TOOL_BUTTON(_radioitem) );
    gtk_tool_button_set_label( GTK_TOOL_BUTTON(toolitem), "Zoom out" );
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Zoom out" );
    gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "toggled",
		      G_CALLBACK(menuitem_tool_change_signal),
		      (gpointer)this );

    // Creating "Zoom fit" button
    pixbuf = gdk_pixbuf_new_from_inline( -1, icon_zoom_fit_inline, FALSE, NULL );
    icon = gtk_image_new_from_pixbuf( pixbuf );
    //toolitem = gtk_radio_tool_button_new_from_widget( GTK_RADIO_TOOL_BUTTON(radioitem) );
    toolitem = gtk_tool_button_new( icon, "Zoom fit" );
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Zoom fit" );
    //gtk_tool_button_set_label( GTK_TOOL_BUTTON(toolitem), "Zoom fit" );
    //gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "clicked",
		      G_CALLBACK(menuitem_zoom_fit_signal),
		      (gpointer)this );
    
    // Creating "Move" button
    pixbuf = gdk_pixbuf_new_from_inline( -1, icon_move_inline, FALSE, NULL );
    icon = gtk_image_new_from_pixbuf( pixbuf );
    toolitem = gtk_radio_tool_button_new_from_widget( GTK_RADIO_TOOL_BUTTON(_radioitem) );
    gtk_tool_button_set_label( GTK_TOOL_BUTTON(toolitem), "Move" );
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Move" );
    gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
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
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Track plot" );
    gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "toggled",
		      G_CALLBACK(menuitem_tool_change_signal),
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
    g_signal_connect( G_OBJECT(_darea), "configure_event",
		      G_CALLBACK(darea_configure_signal), 
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "expose_event",
		      G_CALLBACK(darea_expose_signal), 
		      (gpointer)this );
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


GTKWindow::~GTKWindow()
{

}


void GTKWindow::show( void )
{
    gtk_widget_show_all( _window );
    //
}


void GTKWindow::configure( void )
{
    //std::cout << "Configure\n";
    _width = _darea->allocation.width;
    _height = _darea->allocation.height;

    if( _cairo ) {
	cairo_destroy( _cairo );
	cairo_surface_destroy( _surface );
    }

    _surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, _width, _height );
    _cairo = cairo_create( _surface );

    // Draw
    _frame.set_geometry( _width, _height, 0, 0 );
    _frame.draw( _cairo );
}


gboolean GTKWindow::darea_configure_signal( GtkWidget *widget, 
					     GdkEventConfigure *event, 
					     gpointer object )
{
    GTKWindow *plotter = (GTKWindow *)object;
    plotter->configure();
    return( FALSE );
}


void GTKWindow::expose( int x, int y, int width, int height )
{
    //std::cout << "Expose\n";
    cairo_t *cairo;
    cairo = gdk_cairo_create( _darea->window );
    cairo_set_source_surface( cairo, _surface, 0, 0 );
    cairo_pattern_set_filter( cairo_get_source(cairo), CAIRO_FILTER_FAST );
    cairo_rectangle( cairo, x, y, width, height );
    cairo_clip( cairo );
    cairo_scale( cairo, 50, 50 );
    cairo_translate( cairo, 0, 0 );
    cairo_paint( cairo );
    cairo_destroy( cairo );
}


gboolean GTKWindow::darea_expose_signal( GtkWidget *widget, 
					  GdkEventExpose *event, 
					  gpointer object )
{
    GTKWindow *plotter = (GTKWindow *)object;
    plotter->expose( event->area.x, event->area.y, 
		     event->area.width, event->area.height );
    return( FALSE );
}


void *GTKWindow::build_preferences( GtkWidget *notebook )
{
    // No default action
    return( NULL );
}


void GTKWindow::read_preferences( GtkWidget *notebook, void *pdata )
{
    // No default action
}


void GTKWindow::zoom_fit( void )
{
    //std::cout << "Zoom fit\n";
    double min = -std::numeric_limits<double>::infinity();
    double max = std::numeric_limits<double>::infinity();
    _frame.set_ranges( PLOT_AXIS_X1, min, max );
    _frame.set_ranges( PLOT_AXIS_Y1, min, max );
    //_frame.ruler_autorange_enable( PLOT_AXIS_X1, false, false );
    //_frame.ruler_autorange_enable( PLOT_AXIS_Y1, false, false );
    
    //Redraw and enforce expose
    _frame.draw( _cairo );
    expose( 0, 0, _width, _height );
}


std::string GTKWindow::track_text( double x, double y )
{
    std::string s;
    return( s );
}


// action: 0 for leave, 1 for motion, 2 for enter, 3 for create, 4 for delete
void GTKWindow::track( int action, double x, double y )
{
    double margin[4];
    _frame.get_margins( margin );

    if( action == 4 ) {
	// Delete track window
	if( _trackwindow ) {
	    gtk_widget_destroy( _trackwindow );
	    _trackwindow = NULL;
	}
    } else if( action == 3 ) {
	// Create track window
        _trackwindow = gtk_window_new( GTK_WINDOW_TOPLEVEL );
        gtk_window_set_title( GTK_WINDOW(_trackwindow), "Track" );
        gtk_window_set_policy( GTK_WINDOW(_trackwindow), TRUE, TRUE, FALSE );

        _tracklabel = gtk_label_new( "" );
        gtk_misc_set_alignment( GTK_MISC(_tracklabel), 0.0, 0.0 );
        gtk_label_set_justify( GTK_LABEL(_tracklabel), GTK_JUSTIFY_LEFT );
        gtk_widget_set_size_request( _tracklabel, 300, 150 );
        
        gtk_container_add( GTK_CONTAINER(_trackwindow), _tracklabel );
        g_signal_connect( G_OBJECT(_trackwindow), "delete_event",
                          G_CALLBACK(gtk_true),  // Prevent delete
                          (gpointer)0 );
        gtk_widget_show_all( _trackwindow );
    }

    if( action == 1 || action == 0 ) {
	// Erase old track crosshair
	expose( _end[0], 0, 1, _height );
	expose( 0, _end[1], _width, 1 );
    }

    if( action == 2 || action == 1 ) {
	// Draw new track crosshair
	GdkPoint  p[2];
	GdkColor  color;
	GdkGC    *gc;
	_end[0] = (int)floor(x+0.5);
	_end[1] = (int)floor(y+0.5);
	color.red = color.green = color.blue = 0;
	gc = gdk_gc_new( _darea->window );
	gdk_gc_set_rgb_fg_color( gc, &color );
	if( y >= margin[3] && y <= _height-margin[1] ) {
	    p[0].x = (int)floor(margin[0]+0.5);
	    p[0].y = _end[1];
	    p[1].x = _width-(int)floor(margin[2]+0.5);
	    p[1].y = _end[1];
	    gdk_draw_lines( _darea->window, gc, p, 2 );
	}
	if( x >= margin[0] && x <= _width-margin[2] ) {
	    p[0].x = _end[0];
	    p[0].y = (int)floor(margin[3]+0.5);
	    p[1].x = _end[0];
	    p[1].y = _height-(int)floor(margin[1]+0.5);
	    gdk_draw_lines( _darea->window, gc, p, 2 );
	}
	
	Coordmapper cm = _frame.get_coordmapper( PLOT_AXIS_X1, PLOT_AXIS_Y1 );
	cm.inv_transform( x, y );
	std::string s = track_text( x, y );
	gtk_label_set_text( GTK_LABEL(_tracklabel), s.c_str() );
    }    
}


void GTKWindow::move( int action, double x, double y )
{
    if( action == 0 ) {
	//std::cout << "Move start\n";
	_start[0] = (int)floor(x+0.5);
	_start[1] = (int)floor(y+0.5);
    } else {
	//std::cout << "Move end\n";
	double corners[4];
	int dx[2] = {(int)floor(x+0.5) - _start[0], 
		     (int)floor(y+0.5) - _start[1]};
	double range[4];
	_frame.get_ranges( PLOT_AXIS_X1, range[0], range[2] );
	_frame.get_ranges( PLOT_AXIS_Y1, range[1], range[3] );
	Coordmapper cm = _frame.get_coordmapper( PLOT_AXIS_X1, PLOT_AXIS_Y1 );
	cm.transform( &corners[0], &range[0] );
	cm.transform( &corners[2], &range[2] );
	corners[0] -= dx[0];
	corners[2] -= dx[0];
	corners[1] -= dx[1];
	corners[3] -= dx[1];
	_start[0] += dx[0];
	_start[1] += dx[1];
	cm.inv_transform( &range[0], &corners[0] );
	cm.inv_transform( &range[2], &corners[2] );
	_frame.set_ranges( PLOT_AXIS_X1, range[0], range[2] );
	_frame.set_ranges( PLOT_AXIS_Y1, range[1], range[3] );

	//Redraw and enforce expose
	_frame.draw( _cairo );
	expose( 0, 0, _width, _height );
    }

}


void GTKWindow::zoom_out( double x, double y )
{
    //std::cout << "Zoom out\n";

    double edge[4];
    _frame.get_frame_edges( edge );
    double plotw = edge[2] - edge[0];
    double ploth = edge[1] - edge[3];
    double offsx = (x-edge[0]) / plotw;
    double offsy = (y-edge[3]) / ploth;
    double corners[4];
    corners[0] = x - 2.0*plotw*offsx;
    corners[2] = x + 2.0*plotw*(1.0-offsx);
    corners[1] = y - 2.0*ploth*offsy;
    corners[3] = y + 2.0*ploth*(1.0-offsy);

    // X1 and Y1 axes
    double range[4];
    Coordmapper cm = _frame.get_coordmapper( PLOT_AXIS_X1, PLOT_AXIS_Y1 );
    cm.inv_transform( &range[0], &corners[0] );
    cm.inv_transform( &range[2], &corners[2] );
    _frame.set_ranges( PLOT_AXIS_X1, range[0], range[2] );
    _frame.set_ranges( PLOT_AXIS_Y1, range[1], range[3] );

    // X2 and Y2 axes
    cm = _frame.get_coordmapper( PLOT_AXIS_X2, PLOT_AXIS_Y2 );
    cm.inv_transform( &range[0], &corners[0] );
    cm.inv_transform( &range[2], &corners[2] );
    _frame.set_ranges( PLOT_AXIS_X2, range[0], range[2] );
    _frame.set_ranges( PLOT_AXIS_Y2, range[1], range[3] );

    //Redraw and enforce expose
    _frame.draw( _cairo );
    expose( 0, 0, _width, _height );
}


void GTKWindow::zoom_in( double x, double y )
{
    //std::cout << "Zoom in\n";
    double edge[4];
    _frame.get_frame_edges( edge );
    double plotw = edge[2] - edge[0];
    double ploth = edge[1] - edge[3];
    //std::cout << "plotw = " << plotw << "\n";
    //std::cout << "ploth = " << ploth << "\n";
    double offsx = (x-edge[0]) / plotw;
    double offsy = (y-edge[3]) / ploth;
    //std::cout << "offsx = " << offsx << "\n";
    //std::cout << "offsy = " << offsy << "\n";
    double corners[4];
    corners[0] = x - 0.5*plotw*offsx;
    corners[2] = x + 0.5*plotw*(1.0-offsx);
    corners[1] = y - 0.5*ploth*offsy;
    corners[3] = y + 0.5*ploth*(1.0-offsy);
    //std::cout << "corners[0] = " << corners[0] << "\n";
    //std::cout << "corners[1] = " << corners[1] << "\n";
    //std::cout << "corners[2] = " << corners[2] << "\n";
    //std::cout << "corners[3] = " << corners[3] << "\n";
	
    // X1 and Y1 axes
    double range[4];
    Coordmapper cm = _frame.get_coordmapper( PLOT_AXIS_X1, PLOT_AXIS_Y1 );
    cm.inv_transform( &range[0], &corners[0] );
    cm.inv_transform( &range[2], &corners[2] );
    _frame.set_ranges( PLOT_AXIS_X1, range[0], range[2] );
    _frame.set_ranges( PLOT_AXIS_Y1, range[1], range[3] );
    //std::cout << "range[0] = " << range[0] << "\n";
    //std::cout << "range[1] = " << range[1] << "\n";
    //std::cout << "range[2] = " << range[2] << "\n";
    //std::cout << "range[3] = " << range[3] << "\n";

    // X2 and Y2 axes
    cm = _frame.get_coordmapper( PLOT_AXIS_X2, PLOT_AXIS_Y2 );
    cm.inv_transform( &range[0], &corners[0] );
    cm.inv_transform( &range[2], &corners[2] );
    _frame.set_ranges( PLOT_AXIS_X2, range[0], range[2] );
    _frame.set_ranges( PLOT_AXIS_Y2, range[1], range[3] );
    //std::cout << "range[0] = " << range[0] << "\n";
    //std::cout << "range[1] = " << range[1] << "\n";
    //std::cout << "range[2] = " << range[2] << "\n";
    //std::cout << "range[3] = " << range[3] << "\n";

    //Redraw and enforce expose
    _frame.draw( _cairo );
    expose( 0, 0, _width, _height );
}


void GTKWindow::zoom_window( int action, double x, double y )
{
    if( action == 0 ) {
	//std::cout << "Zoom window start\n";
	_start[0] = _end[0] = (int)floor(x+0.5);
	_start[1] = _end[1] = (int)floor(y+0.5);
    } else if( action == 1 ) {
	//std::cout << "Zoom window modify\n";
	// Erase old zoom box
	int x0, y0, width, height;
	x0 = _start[0] < _end[0] ? _start[0] : _end[0];
	y0 = _start[1] < _end[1] ? _start[1] : _end[1];
	width = abs( _start[0] - _end[0] );
	height = abs( _start[1] - _end[1] );
	expose( x0, y0, 1, height+1 );
	expose( x0, y0, width+1, 1 );
	expose( x0+width, y0, 1, height+1 );
	expose( x0, y0+height, width+1, 1 );

	// Draw new zoom box
	GdkPoint  p[5];
	GdkColor  color;
	GdkGC    *gc;
	_end[0] = (int)floor(x+0.5);
	_end[1] = (int)floor(y+0.5);
	color.red = color.green = color.blue = 0;
	gc = gdk_gc_new( _darea->window );
	gdk_gc_set_rgb_fg_color( gc, &color );
	p[0].x = _start[0];
	p[0].y = _start[1];
	p[1].x = _end[0];
	p[1].y = _start[1];
	p[2].x = _end[0];
	p[2].y = _end[1];
	p[3].x = _start[0];
	p[3].y = _end[1];
	p[4].x = _start[0];
	p[4].y = _start[1];
	gdk_draw_lines( _darea->window, gc, p, 5 );

    } else {
	//std::cout << "Zoom window end\n";
	double corners[4];
	double range[4];
	corners[0] = _start[0] < _end[0] ? _start[0] : _end[0];
	corners[2] = _start[0] > _end[0] ? _start[0] : _end[0];
	corners[1] = _start[1] > _end[1] ? _start[1] : _end[1];
	corners[3] = _start[1] < _end[1] ? _start[1] : _end[1];

	// X1 and Y1 axes
	Coordmapper cm = _frame.get_coordmapper( PLOT_AXIS_X1, PLOT_AXIS_Y1 );
	cm.inv_transform( &range[0], &corners[0] );
	cm.inv_transform( &range[2], &corners[2] );
	_frame.set_ranges( PLOT_AXIS_X1, range[0], range[2] );
	_frame.set_ranges( PLOT_AXIS_Y1, range[1], range[3] );

	// X2 and Y2 axes
	cm = _frame.get_coordmapper( PLOT_AXIS_X2, PLOT_AXIS_Y2 );
	cm.inv_transform( &range[0], &corners[0] );
	cm.inv_transform( &range[2], &corners[2] );
	_frame.set_ranges( PLOT_AXIS_X2, range[0], range[2] );
	_frame.set_ranges( PLOT_AXIS_Y2, range[1], range[3] );

	//Redraw and enforce expose
	_frame.draw( _cairo );
	expose( 0, 0, _width, _height );
    }
}


void GTKWindow::darea_motion( GdkEventMotion *event ) 
{
    if( (_tool == TOOL_MOVE && (event->state & GDK_BUTTON1_MASK)) || event->state & GDK_BUTTON2_MASK ) {
	move( 1, event->x, event->y );
    } else if( _tool == TOOL_TRACK ) {
	track( 1, event->x, event->y );
    } else if( _tool == TOOL_ZOOM_IN && (event->state & GDK_BUTTON1_MASK) ) {
	zoom_window( 1, event->x, event->y );
    } else if( event->state & GDK_BUTTON2_MASK ) {
	zoom_window( 1, event->x, event->y );
    }

    // Statusbar
    Coordmapper cm = _frame.get_coordmapper( PLOT_AXIS_X1, PLOT_AXIS_Y1 );
    double cx = event->x;
    double cy = event->y;
    cm.inv_transform( cx, cy );
    std::stringstream ss;
    ss << "("<< cx << ", " << cy << ")";
    gtk_statusbar_pop( GTK_STATUSBAR(_statusbar), 0 );
    gtk_statusbar_push( GTK_STATUSBAR(_statusbar), 0, ss.str().c_str() );

    int x, y;
    GdkModifierType state;
    gdk_window_get_pointer( event->window, &x, &y, &state );
}


void GTKWindow::darea_enter( GdkEventCrossing *event )
{
    if( _tool == TOOL_TRACK )
	track( 2, 0, 0 );
}


void GTKWindow::darea_leave( GdkEventCrossing *event )
{
    if( _tool == TOOL_TRACK )
	track( 0, event->x, event->y );

    // Statusbar
    gtk_statusbar_pop( GTK_STATUSBAR(_statusbar), 0 );
    gtk_statusbar_push( GTK_STATUSBAR(_statusbar), 0, "" );
}


gboolean GTKWindow::darea_motion_signal( GtkWidget *widget, 
					 GdkEventMotion *event,
					 gpointer object )
{
    GTKWindow *plotter = (GTKWindow *)object;
    plotter->darea_motion( event );
    return( FALSE );  
}


gboolean GTKWindow::darea_enter_signal( GtkWidget *widget, 
					GdkEventCrossing *event,
					gpointer object )
{
    GTKWindow *plotter = (GTKWindow *)object;
    plotter->darea_enter( event );
    return( FALSE );  
}


gboolean GTKWindow::darea_leave_signal( GtkWidget *widget, 
					GdkEventCrossing *event,
					gpointer object )
{
    GTKWindow *plotter = (GTKWindow *)object;
    plotter->darea_leave( event );
    return( FALSE );
}


void GTKWindow::darea_button( GdkEventButton *event )
{
    if( (_tool == TOOL_MOVE && event->type == GDK_BUTTON_PRESS && event->button == 1) ||
	(event->type == GDK_BUTTON_PRESS && event->button == 2) ) {
	move( 0, event->x, event->y );
    } else if( (_tool == TOOL_MOVE && event->type == GDK_BUTTON_RELEASE && event->button == 1) || 
	       (event->type == GDK_BUTTON_RELEASE && event->button == 2) ) {
	move( 2, event->x, event->y );
    } else if( _tool == TOOL_ZOOM_OUT && event->type == GDK_BUTTON_RELEASE && event->button == 1 ) {
	zoom_out( event->x, event->y );
    } else if( _tool == TOOL_ZOOM_IN && event->type == GDK_BUTTON_PRESS && event->button == 1 ) {
	zoom_window( 0, event->x, event->y );
    } else if( _tool == TOOL_ZOOM_IN && event->type == GDK_BUTTON_RELEASE && event->button == 1 ) {
	if( fabs( _start[0] - event->x ) < 2.0 ||
	    fabs( _start[1] - event->y ) < 2.0 )
	    zoom_in( event->x, event->y );
	else
	    zoom_window( 2, event->x, event->y );
    }
}
					     

gboolean GTKWindow::darea_button_signal( GtkWidget *widget, 
					 GdkEventButton *event,
					 gpointer object )
{
    GTKWindow *plotter = (GTKWindow *)object;
    plotter->darea_button( event );
    return( FALSE );
}


void GTKWindow::delete_window( void )
{
    if( _tool == TOOL_TRACK )
	track( 4, 0, 0 );

    _plotter->delete_window( this );
}


gboolean GTKWindow::window_delete_signal( GtkWidget *widget, 
					  GdkEventExpose *event, 
					  gpointer object )
{
    GTKWindow *plotter = (GTKWindow *)object;
    plotter->delete_window();
    return( FALSE );
}


void GTKWindow::menuitem_preferences( GtkMenuItem *menuitem )
{
    std::cout << "menu preferences\n";
    GTKPreferences preferences( this, _window, &_frame );
    preferences.run();
}


void GTKWindow::menuitem_preferences_signal( GtkMenuItem *menuitem,
					     gpointer object )
{
    GTKWindow *plotter = (GTKWindow *)object;
    plotter->menuitem_preferences( menuitem );
}


void GTKWindow::menuitem_quit_signal( GtkMenuItem *menuitem,
				      gpointer object )
{
    GTKWindow *plotter = (GTKWindow *)object;
    plotter->delete_window();
}


void GTKWindow::menuitem_tool_change( GtkToolButton *button )
{
    int tool;
    const char *label = gtk_tool_button_get_label( button );
    if( !strcmp( label, "Zoom in" ) ) 
	tool = TOOL_ZOOM_IN;
    else if( !strcmp( label, "Zoom out" ) ) 
	tool = TOOL_ZOOM_OUT;
    else if( !strcmp( label, "Move" ) ) 
	tool = TOOL_MOVE;
    else if( !strcmp( label, "Track" ) )
	tool = TOOL_TRACK;
    else {
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


void GTKWindow::menuitem_tool_change_signal( GtkToolButton *button,
					     gpointer object )
{
    GTKWindow *plotter = (GTKWindow *)object;
    plotter->menuitem_tool_change( button );
}


void GTKWindow::menuitem_hardcopy_signal( GtkToolButton *button,
					  gpointer object )
{
    GTKWindow *plotter = (GTKWindow *)object;
    plotter->hardcopy();
}


void GTKWindow::menuitem_zoom_fit_signal( GtkToolButton *button,
					  gpointer object )
{
    GTKWindow *plotter = (GTKWindow *)object;
    plotter->zoom_fit();
}


void GTKWindow::draw_and_expose( void )
{
    //Redraw and enforce expose
    _frame.draw( _cairo );
    expose( 0, 0, _width, _height );
}


void GTKWindow::hardcopy( void )
{
    GTKHardcopy hardcopy( _window, &_frame, _width, _height );
    hardcopy.run();

    // Reset frame settings
    _frame.set_geometry( _width, _height, 0, 0 );
}

















