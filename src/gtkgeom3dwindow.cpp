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


#include <gtk/gtkgl.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "gtkgeom3dwindow.hpp"
#include "camera.hpp"
#include "icons.hpp"


#define TOOL_UNKNOWN  -1
#define TOOL_MOVE      0
#define TOOL_ZOOM_IN   1
#define TOOL_ZOOM_OUT  2
#define TOOL_TRACK     3


GTKGeom3DWindow::GTKGeom3DWindow( GTKPlotter &plotter,
				  const Geometry &geom )
    : _plotter(plotter), _geom(geom), _width(640), _height(480), _camera(NULL)
{
    init_camera_and_model();
    init_window();

    // Initialize OpenGL
    std::cout << "Init OpenGL config\n";
    GdkGLConfigMode mode = (GdkGLConfigMode)( GDK_GL_MODE_RGBA |
					      GDK_GL_MODE_DEPTH |
					      GDK_GL_MODE_DOUBLE );
    GdkGLConfig *gl_config = gdk_gl_config_new_by_mode( mode );
    if( !gl_config )
        g_assert_not_reached();

    if( !gtk_widget_set_gl_capability( _darea, gl_config, NULL, TRUE,
				       GDK_GL_RGBA_TYPE ) )
        g_assert_not_reached();

    // Must be set after set_gl_capability
    g_signal_connect( G_OBJECT(_darea), "configure_event",
		      G_CALLBACK(darea_configure_signal), 
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "expose_event",
		      G_CALLBACK(darea_expose_signal), 
		      (gpointer)this );

    gtk_widget_show_all( _window );
}


void GTKGeom3DWindow::init_camera_and_model( void )
{
    // Calculate scale
    _center = 0.5*(_geom.origo() + _geom.max());
    Vec3D size = _geom.max() - _geom.origo();
    _scale = 1.0/size.norm2();

    // Init camera
    _camera = new PerspectiveCamera;
    _camera->set_zplanes( 1.0, 3.0 );
    _camera->set_camera_location( Vec3D(0,0,2) );
    _camera->set_target_location( Vec3D(0,0,0) );
    _camera->set_camera_up( Vec3D(0,1,0) );
    _camera->set_size( _width, _height );
    _camera->set_field_of_view( 0.25 );

    // Init model transformation
    _modeltrans = Transformation::unity();
}


void GTKGeom3DWindow::init_window( void )
{
    std::cout << "Init Window\n";

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
    /*
    g_signal_connect( G_OBJECT(item_hardcopy), "activate",
		      G_CALLBACK(menuitem_hardcopy_signal),
		      (gpointer)this );
    */
    gtk_box_pack_start( GTK_BOX(vbox), _menubar, FALSE, TRUE, 0 );

    // Add Edit/Configure menu
    GtkWidget *menu_edit, *item_edit, *item_preferences;
    menu_edit = gtk_menu_new();
    item_preferences = gtk_menu_item_new_with_mnemonic( "_Preferences" );
    gtk_menu_shell_append( GTK_MENU_SHELL(menu_edit), item_preferences );
    item_edit = gtk_menu_item_new_with_mnemonic( "_Edit" );
    gtk_menu_item_set_submenu( GTK_MENU_ITEM(item_edit), menu_edit );
    gtk_menu_shell_insert( GTK_MENU_SHELL(_menubar), item_edit, 1 );
    /*
    g_signal_connect( G_OBJECT(item_preferences), "activate",
                      G_CALLBACK(menuitem_preferences_signal),
                      (gpointer)this );
    */

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
    /*
    g_signal_connect( G_OBJECT(toolitem), "clicked",
		      G_CALLBACK(menuitem_hardcopy_signal),
		      (gpointer)this );
    */
    
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
    /*
    g_signal_connect( G_OBJECT(_darea), "configure_event",
		      G_CALLBACK(darea_configure_signal), 
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "expose_event",
		      G_CALLBACK(darea_expose_signal), 
		      (gpointer)this );
    */
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
    if( _camera )
	delete _camera;
}


void GTKGeom3DWindow::delete_window( void )
{
    _plotter.delete_window( this );
}


void GTKGeom3DWindow::configure( void )
{
    std::cout << "Configure\n";

    // Initialize OpenGL context
    GdkGLContext *gl_context = gtk_widget_get_gl_context( _darea );
    GdkGLDrawable *gl_drawable = gtk_widget_get_gl_drawable( _darea );
    if( !gdk_gl_drawable_gl_begin( gl_drawable, gl_context ) )
	g_assert_not_reached();

    // Set OpenGL viewport
    GtkAllocation alloc;
    gtk_widget_get_allocation( _darea, &alloc );
    glViewport( 0, 0, alloc.width, alloc.height );
    _width = alloc.width;
    _height = alloc.height;

    // Set camera
    _camera->set_size( _width, _height );

    // Setup OpenGL rendering
    float light_ambient[] = { 1.0, 1.0, 1.0, 0.0 };
    float light_position[] = { 0.0, 0.0, -50.0, 0.0 };
    glLightfv( GL_LIGHT0, GL_AMBIENT, light_ambient );
    glLightfv( GL_LIGHT0, GL_POSITION, light_position );

    glEnable( GL_LIGHTING );
    glEnable( GL_LIGHT0 );
    glEnable( GL_DEPTH_TEST );

    // Close OpenGL context
    gdk_gl_drawable_gl_end( gl_drawable );
}


void GTKGeom3DWindow::draw_bbox( void )
{
    Vec3D min = _scale*(_geom.origo()-_center);
    Vec3D max = _scale*(_geom.max()-_center);

    glDisable( GL_LIGHTING );
    glColor3f( 0.0, 0.0, 0.0 );

    glBegin( GL_LINE_LOOP );
    glVertex3d( min[0], min[1], min[2] );
    glVertex3d( min[0], max[1], min[2] );
    glVertex3d( max[0], max[1], min[2] );
    glVertex3d( max[0], min[1], min[2] );
    glEnd();

    glBegin( GL_LINE_LOOP );
    glVertex3d( min[0], min[1], max[2] );
    glVertex3d( min[0], max[1], max[2] );
    glVertex3d( max[0], max[1], max[2] );
    glVertex3d( max[0], min[1], max[2] );
    glEnd();

    glBegin( GL_LINES );
    glVertex3d( min[0], min[1], min[2] );
    glVertex3d( min[0], min[1], max[2] );

    glVertex3d( min[0], max[1], min[2] );
    glVertex3d( min[0], max[1], max[2] );

    glVertex3d( max[0], max[1], min[2] );
    glVertex3d( max[0], max[1], max[2] );

    glVertex3d( max[0], min[1], min[2] );
    glVertex3d( max[0], min[1], max[2] );
    glEnd();
}


void GTKGeom3DWindow::draw_model( void )
{
    glEnable( GL_LIGHTING );
    float material_diffuse[] =  { 0.2, 0.2, 0.8, 0.0 };
    float material_ambient[] =  { 0.0, 0.0, 0.2, 0.0 };
    glMaterialfv( GL_FRONT, GL_DIFFUSE, material_diffuse );
    glMaterialfv( GL_FRONT, GL_AMBIENT, material_ambient );

    Vec3D x[3];
    for( int32_t i = 0; i < _geom.surface_trianglec(); i++ ) {
	const VTriangle &tri = _geom.surface_triangle(i);
	
	uint32_t v0 = tri[0];
	x[0] = _scale*(_geom.surface_vertex(v0)-_center);
	uint32_t v1 = tri[1];
	x[1] = _scale*(_geom.surface_vertex(v1)-_center);
	uint32_t v2 = tri[2];
	x[2] = _scale*(_geom.surface_vertex(v2)-_center);

	Vec3D norm = cross( x[1]-x[0], x[2]-x[0] );
	norm.normalize();

	glBegin( GL_TRIANGLES );
	glNormal3dv( &norm[0] );
	glVertex3dv( &x[0][0] );
	glVertex3dv( &x[1][0] );
	glVertex3dv( &x[2][0] );
	glEnd();
    }
}


void GTKGeom3DWindow::draw_cut_planes( void )
{
    glEnable( GL_LIGHTING );
    float material_diffuse[] =  { 0.2, 0.2, 0.8, 0.0 };
    float material_ambient[] =  { 0.0, 0.0, 0.2, 0.0 };
    glMaterialfv( GL_FRONT, GL_DIFFUSE, material_diffuse );
    glMaterialfv( GL_FRONT, GL_AMBIENT, material_ambient );

    for( int32_t p = 0; p < 6; p++ ) {
	for( uint32_t i = 0; i < _csurface[p].size(); i+=3 ) {
	
	    Vec3D x0 = _scale*(_csurface[p][i+0]-_center);
	    Vec3D x1 = _scale*(_csurface[p][i+1]-_center);
	    Vec3D x2 = _scale*(_csurface[p][i+2]-_center);

	    Vec3D norm = cross( x1-x0, x2-x0 );
	    norm.normalize();
	
	    glBegin( GL_TRIANGLES );
	    glNormal3dv( &norm[0] );
	    glVertex3dv( &x0[0] );
	    glVertex3dv( &x1[0] );
	    glVertex3dv( &x2[0] );
	    glEnd();
	}
    }
}


int32_t GTKGeom3DWindow::case2d( const int i[3], const int vb[3] )
{
    int res = 0;
    uint32_t node;
    int32_t j[3] = { i[0], i[1], i[2] };

    // Node 1 (x,y)
    node = _geom.mesh((j[2]*_geom.size(1) + j[1])*_geom.size(0) + j[0]);
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += 1;

    // Node 2 (x+1,y)
    j[vb[0]]++;
    node = _geom.mesh((j[2]*_geom.size(1) + j[1])*_geom.size(0) + j[0]);
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += 2;

    // Node 3 (x+1,y+1)
    j[vb[1]]++;
    node = _geom.mesh((j[2]*_geom.size(1) + j[1])*_geom.size(0) + j[0]);
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += 4;

    // Node 4 (x,y+1)
    j[vb[0]]--;
    node = _geom.mesh((j[2]*_geom.size(1) + j[1])*_geom.size(0) + j[0]);
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += 8;

    return( res );
}


void GTKGeom3DWindow::cplane_add_vertex( int32_t p, const int32_t i[3], const int32_t vb[3], double dx, double dy )
{
    double u[3] = { i[0], i[1], i[2] };
    u[vb[0]] += dx;
    u[vb[1]] += dy;

    Vec3D x( _geom.origo(0)+_geom.h()*u[0], _geom.origo(1)+_geom.h()*u[1], _geom.origo(2)+_geom.h()*u[2] );
    _csurface[p].push_back( x );
}


double GTKGeom3DWindow::cplane_dist( const int32_t i[3], const int32_t vb[3], int32_t dx, int32_t dy, int32_t dir )
{
    int32_t j[3] = { i[0], i[1], i[2] };
    j[vb[0]] += dx;
    j[vb[1]] += dy;
    return( _geom.solid_dist( j[0], j[1], j[2], dir )/255.0 );
}


void GTKGeom3DWindow::build_cut_plane( int32_t p, int32_t const vb[3], int32_t level )
{
    int32_t i[3];
    i[vb[2]] = level;

    // Go through mesh squares in cut plane
    int32_t sizex = _geom.size(vb[0])-1;
    int32_t sizey = _geom.size(vb[1])-1;
    for( int32_t y = 0; y < sizey; y++ ) {
	i[vb[1]] = y;

	for( int32_t x = 0; x < sizex; x++ ) {
	    i[vb[0]] = x;

	    // Construct case number
	    int32_t cn = case2d( i, vb );

	    double dist;
	    double dist2;
	    switch( cn ) {
	    case 0:
		// Fully outside
		break;
	    case 1:
		// (i,j) in
		dist = cplane_dist( i, vb, 0, 1, 2*vb[0] );
		cplane_add_vertex( p, i, vb, dist, 0.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		break;
	    case 2:
		// (i+1,j) in
		dist = cplane_dist( i, vb, 0, 0, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist, 0.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		break;
	    case 3:
		// (i,j) and (i+1,j) in
		dist = cplane_dist( i, vb, 0, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 1.0, 1.0-dist );
		break;
	    case 4:
		// (i+1,j+1) in
		dist = cplane_dist( i, vb, 0, 1, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist, 1.0 );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 1.0, dist );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		break;
	    case 5:
		// (i,j) and (i+1,j+1) in
		dist = cplane_dist( i, vb, 0, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 0.0 );
		dist = cplane_dist( i, vb, 0, 1, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist, 1.0 );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 1.0, dist );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		break;
	    case 6:
		// (i+1,j) and (i+1,j+1) in
		dist = cplane_dist( i, vb, 0, 1, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist, 1.0 );
		dist = cplane_dist( i, vb, 0, 0, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, dist, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		break;
	    case 7:
		// (i,j), (i+1,j) and (i+1,j+1) in
		dist = cplane_dist( i, vb, 0, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		dist2 = cplane_dist( i, vb, 0, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, dist2, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, dist2, 1.0 );
		break;
	    case 8:
		dist = cplane_dist( i, vb, 0, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 0.0, dist );
		dist = cplane_dist( i, vb, 1, 1, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		break;
	    case 9:
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 1, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 1.0 );
		cplane_add_vertex( p, i, vb, 1.0-dist, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 0.0 );
		break;
	    case 10:
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		dist = cplane_dist( i, vb, 0, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 0.0, dist );
		dist = cplane_dist( i, vb, 1, 1, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 1.0 );
		dist = cplane_dist( i, vb, 0, 0, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 1.0, 1.0-dist );
		break;
	    case 11:
		dist = cplane_dist( i, vb, 1, 1, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		dist2 = cplane_dist( i, vb, 1, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 1.0, 1.0-dist2 );
		cplane_add_vertex( p, i, vb, 1.0-dist, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0-dist2 );
		break;
	    case 12:
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		dist = cplane_dist( i, vb, 0, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 0.0, dist );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 1.0, dist );
		cplane_add_vertex( p, i, vb, 1.0, dist );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		break;
	    case 13:
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0-dist, 0.0 );
		dist2 = cplane_dist( i, vb, 1, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 1.0, dist2 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 1.0, dist2 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		break;
	    case 14:
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		dist = cplane_dist( i, vb, 0, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 0.0, dist );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, dist );
		dist2 = cplane_dist( i, vb, 0, 0, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist2, 0.0 );
		cplane_add_vertex( p, i, vb, dist2, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		break;
	    case 15:
		// Fully inside
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		break;
	    }
	    
	}
    }
}


void GTKGeom3DWindow::build_cut_planes( void )
{
    int vb[3];

    // Z=0
    vb[0] = 0;
    vb[1] = 1;
    vb[2] = 2;
    build_cut_plane( 0, vb, 0 );

    // Z=size(2)-1
    vb[0] = 1;
    vb[1] = 0;
    vb[2] = 2;
    build_cut_plane( 1, vb, _geom.size(2)-1 );
}


void GTKGeom3DWindow::draw( void )
{
    // Clear
    glClearColor( 1.0, 1.0, 1.0, 0.0 );
    glClearDepth( 1.0 );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    // Set model transformation
    glMultMatrixd( &_modeltrans[0] );

    // Draw model
    build_cut_planes();
    draw_cut_planes();
    draw_model();

    // Draw bbox
    draw_bbox();
}


void GTKGeom3DWindow::expose( void )
{
    // Initialize OpenGL context
    GdkGLContext *gl_context = gtk_widget_get_gl_context( _darea );
    GdkGLDrawable *gl_drawable = gtk_widget_get_gl_drawable( _darea );
    if( !gdk_gl_drawable_gl_begin( gl_drawable, gl_context ) )
	g_assert_not_reached();

    // Setup camera and draw
    _camera->gl_initalize_camera();
    draw();

    // Finish draw and close OpenGL context
    if( gdk_gl_drawable_is_double_buffered( gl_drawable) )
        gdk_gl_drawable_swap_buffers( gl_drawable ); 
    else
        glFlush();
    gdk_gl_drawable_gl_end( gl_drawable );
}


// action 0: button press (move start), 1: pointer move, 2: release button (move end)
void GTKGeom3DWindow::move( int action, double x, double y )
{
    if( action == 0 ) {
	_oldx = x;
	_oldy = y;
    } else {
	double dx = x-_oldx;
	double dy = y-_oldy;
	_modeltrans.rotate_y_before( -3.0*dx/_width );
	_modeltrans.rotate_x_before( -3.0*dy/_height );
	_oldx = x;
	_oldy = y;
	gtk_widget_queue_draw_area( _darea, 0, 0, _width, _height );
    }
}


void GTKGeom3DWindow::zoom_out( double x, double y )
{
    _camera->set_view_relative( x, y, 1.414 );
    gtk_widget_queue_draw_area( _darea, 0, 0, _width, _height );
}


void GTKGeom3DWindow::zoom_in( double x, double y )
{
    _camera->set_view_relative( x, y, 0.707 );
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
    init_camera_and_model();
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


void GTKGeom3DWindow::menuitem_hardcopy_signal( GtkToolButton *button,
						gpointer object )
{
    //GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
}


void GTKGeom3DWindow::menuitem_zoom_fit_signal( GtkToolButton *button,
						gpointer object )
{
    GTKGeom3DWindow *window = (GTKGeom3DWindow *)object;
    window->zoom_fit();
}

