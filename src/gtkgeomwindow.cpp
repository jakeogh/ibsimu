/*! \file gtkgeomwindow.cpp
 *  \brief Source code for gtkgeomwindow.cpp
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
#include "gtkgeomwindow.hpp"
#include "gtkparticlediagdialog.hpp"
#include "gtkfielddiagdialog.hpp"
#include "epot_efield.hpp"
#include "icons.hpp"


#define TOOL_UNKNOWN  -1
#define TOOL_PARTICLE_DIAG 4
#define TOOL_FIELD_DIAG 5


GTKGeomWindow::GTKGeomWindow( class GTKPlotter       *plotter, 
			      const Geometry         *geom,
			      const ScalarField      *epot,
			      const ScalarField      *scharge,
			      const VectorField      *bfield,
			      const ParticleDataBase *pdb )
    : GTKWindow(plotter), _geomplot(&_frame,geom), 
      _geom(geom), _epot(epot), _scharge(scharge), _bfield(bfield), _pdb(pdb)
{
    std::cout << "GTKGeomWindow constructor\n";

    // Setup GeomPlot
    _geomplot.set_epot( epot );
    _geomplot.set_scharge( scharge );
    _geomplot.set_particle_database( pdb );

    // Set window title
    gtk_window_set_title( GTK_WINDOW(_window), "Simulation geometry" );

    // Adding geometry window specific tools to toolbar
    // Creating separator
    GtkToolItem *toolitem = gtk_separator_tool_item_new();
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );

    // Creating "Particle diagnostics" button
    GdkPixbuf *pixbuf = gdk_pixbuf_new_from_inline( -1, icon_particle_diag_inline, FALSE, NULL );
    GtkWidget *icon = gtk_image_new_from_pixbuf( pixbuf );
    toolitem = gtk_radio_tool_button_new_from_widget( GTK_RADIO_TOOL_BUTTON(_radioitem) );
    gtk_tool_button_set_label( GTK_TOOL_BUTTON(toolitem), "Particle diagnostics" );
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Particle diagnostics" );
    gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "toggled",
		      G_CALLBACK(menuitem_tool_change_signal),
		      (gpointer)this );

    // Creating "Field diagnostics" button
    pixbuf = gdk_pixbuf_new_from_inline( -1, icon_field_diag_inline, FALSE, NULL );
    icon = gtk_image_new_from_pixbuf( pixbuf );
    toolitem = gtk_radio_tool_button_new_from_widget( GTK_RADIO_TOOL_BUTTON(_radioitem) );
    gtk_tool_button_set_label( GTK_TOOL_BUTTON(toolitem), "Field diagnostics" );
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Field diagnostics" );
    gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "toggled",
		      G_CALLBACK(menuitem_tool_change_signal),
		      (gpointer)this );

    // Creating separator
    toolitem = gtk_separator_tool_item_new();
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );

    // Creating view combobox and level spinbutton
    _combobox = gtk_combo_box_new_text();
    gtk_combo_box_append_text( GTK_COMBO_BOX(_combobox), "XY" );
    if( geom->geom_mode() == MODE_3D ) {
	gtk_combo_box_append_text( GTK_COMBO_BOX(_combobox), "XZ" );
        gtk_combo_box_append_text( GTK_COMBO_BOX(_combobox), "YZ" );
	gtk_combo_box_set_active( GTK_COMBO_BOX(_combobox), 0 );
	_spinbutton = gtk_spin_button_new_with_range( 0, geom->size(2)-1, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON(_spinbutton), geom->size(2)/2 );
	_geomplot.set_view( VIEW_XY, geom->size(2)/2 );
    } else {
        gtk_combo_box_set_active( GTK_COMBO_BOX(_combobox), 0 );
        _spinbutton = gtk_spin_button_new_with_range( 0, 0, 1 );
        gtk_spin_button_set_value( GTK_SPIN_BUTTON(_spinbutton), 0 );
    }
    toolitem = gtk_tool_item_new();
    gtk_container_add( GTK_CONTAINER(toolitem), _combobox );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Select view direction" );

    gtk_spin_button_set_digits( GTK_SPIN_BUTTON(_spinbutton), 0 );
    toolitem = gtk_tool_item_new();
    gtk_container_add( GTK_CONTAINER(toolitem), _spinbutton );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Select view level" );
    g_signal_connect( G_OBJECT(_combobox), "changed",
                      G_CALLBACK(combobox_signal),
                      (gpointer)this );
    g_signal_connect( G_OBJECT(_spinbutton), "value-changed",
                      G_CALLBACK(spinbutton_signal),
                      (gpointer)this );

    // Drawing area signals
    g_signal_connect( G_OBJECT(_darea), "button_press_event",
		      G_CALLBACK(darea_button_signal2),
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "button_release_event",
		      G_CALLBACK(darea_button_signal2),
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "motion_notify_event",
		      G_CALLBACK(darea_motion_signal2),
		      (gpointer)this );

    update_view();
    show();
}


void GTKGeomWindow::field_diag( int action, double x, double y )
{
    int x0, y0, width = 0, height = 0;

    if( action == 0 ) {
	// Start
	_start[0] = _end[0] = (int)floor(x+0.5);
	_start[1] = _end[1] = (int)floor(y+0.5);

    } else if( action == 1 || action == 2 ) {
	// Erase old line
	x0 = _start[0] < _end[0] ? _start[0] : _end[0];
	y0 = _start[1] < _end[1] ? _start[1] : _end[1];
	width = abs( _start[0] - _end[0] );
	height = abs( _start[1] - _end[1] );
	expose( x0, y0, width+1, height+1 );
    }

    if( action == 1 ) {
	// Draw new line
	_end[0] = (int)floor(x+0.5);
	_end[1] = (int)floor(y+0.5);
	GdkColor  color;
	GdkGC    *gc;
	color.red = color.green = color.blue = 0;
	gc = gdk_gc_new( _darea->window );
	gdk_gc_set_rgb_fg_color( gc, &color );
	gdk_draw_line( _darea->window, gc, _start[0], _start[1], _end[0], _end[1] );

    } else if( action == 2 ) {
	// Done, start field diagnostics dialog
	Coordmapper cm = _frame.get_coordmapper( PLOT_AXIS_X1, PLOT_AXIS_Y1 );
	double x1[2] = {_start[0], _start[1]};
	double x2[2] = {_end[0], _end[1]};
	cm.inv_transform( x1[0], x1[1] );
	cm.inv_transform( x2[0], x2[1] );
	double p1[3];
	double p2[3];
	p1[_geomplot.vb(0)] = x1[0];
	p2[_geomplot.vb(0)] = x2[0];
	p1[_geomplot.vb(1)] = x1[1];
	p2[_geomplot.vb(1)] = x2[1];
	p1[_geomplot.vb(2)] = _geomplot.get_level_si();
	p2[_geomplot.vb(2)] = _geomplot.get_level_si();
	GTKFieldDiagDialog fielddiag( _window, _plotter, p1, p2 );
	fielddiag.run();
    }
}


void GTKGeomWindow::particle_diag( int action, double x, double y )
{
    int x0, y0, width = 0, height = 0;

    if( action == 0 ) {
	// Start
	_start[0] = _end[0] = (int)floor(x+0.5);
	_start[1] = _end[1] = (int)floor(y+0.5);
    } else if( action == 1 || action == 2 ) {
	// Erase old line
	x0 = _start[0] < _end[0] ? _start[0] : _end[0];
	y0 = _start[1] < _end[1] ? _start[1] : _end[1];
	width = abs( _start[0] - _end[0] );
	height = abs( _start[1] - _end[1] );
	if( width > height )
	    expose( x0, _start[1], width+1, 1 );
	else
	    expose( _start[0], y0, 1, height+1 );
    }

    if( action == 1 ) {
	// Draw new line
	_end[0] = (int)floor(x+0.5);
	_end[1] = (int)floor(y+0.5);
	x0 = _start[0] < _end[0] ? _start[0] : _end[0];
	y0 = _start[1] < _end[1] ? _start[1] : _end[1];
	width = abs( _start[0] - _end[0] );
	height = abs( _start[1] - _end[1] );
	GdkColor  color;
	GdkGC    *gc;
	color.red = color.green = color.blue = 0;
	gc = gdk_gc_new( _darea->window );
	gdk_gc_set_rgb_fg_color( gc, &color );
	if( width > height )
	    gdk_draw_line( _darea->window, gc, _start[0], _start[1], 
			   _end[0], _start[1] );
	else
	    gdk_draw_line( _darea->window, gc, _start[0], _start[1], 
			   _start[0], _end[1] );
    }

    if( action == 2 ) {
	// If no particles
	if( _pdb == NULL ) {
	    GtkWidget *dialog = gtk_message_dialog_new( GTK_WINDOW(_window),
							GTK_DIALOG_DESTROY_WITH_PARENT,
							GTK_MESSAGE_ERROR,
							GTK_BUTTONS_CLOSE,
							"No particle database found." );
	    gtk_dialog_run( GTK_DIALOG(dialog) );
	    gtk_widget_destroy( dialog );
	    return;
	}

	// Done, start particle diagnostics dialog
	int crd;
	double val;
	Coordmapper cm = _frame.get_coordmapper( PLOT_AXIS_X1, PLOT_AXIS_Y1 );
	double c[2] = {_start[0], _start[1]};
	cm.inv_transform( c[0], c[1] );
	if( width > height ) {
	    crd = _geomplot.vb(1);
	    val = c[1];
	} else {
	    crd = _geomplot.vb(0);
	    val = c[0];
	}
	GTKParticleDiagDialog particlediag( _window, _plotter, crd, val );
	particlediag.run();
    }
}


void GTKGeomWindow::darea_motion2( GdkEventMotion *event )
{
    if( _tool == TOOL_PARTICLE_DIAG && (event->state & GDK_BUTTON1_MASK) ) {
	particle_diag( 1, event->x, event->y );
    } else if( _tool == TOOL_FIELD_DIAG && (event->state & GDK_BUTTON1_MASK) ) {
	field_diag( 1, event->x, event->y );
    }
}


void GTKGeomWindow::darea_button2( GdkEventButton *event )
{
    if( _tool == TOOL_PARTICLE_DIAG && event->type == GDK_BUTTON_PRESS && event->button == 1 ) {
	particle_diag( 0, event->x, event->y );
    } else if( _tool == TOOL_PARTICLE_DIAG && event->type == GDK_BUTTON_RELEASE && event->button == 1 ) {
	particle_diag( 2, event->x, event->y );
    } else if( _tool == TOOL_FIELD_DIAG && event->type == GDK_BUTTON_PRESS && event->button == 1 ) {
	field_diag( 0, event->x, event->y );
    } else if( _tool == TOOL_FIELD_DIAG && event->type == GDK_BUTTON_RELEASE && event->button == 1 ) {
	field_diag( 2, event->x, event->y );
    }
}


gboolean GTKGeomWindow::darea_motion_signal2( GtkWidget *widget, 
					      GdkEventMotion *event,
					      gpointer object )
{
    GTKGeomWindow *plotter = (GTKGeomWindow *)object;
    plotter->darea_motion2( event );
    return( FALSE );
}


gboolean GTKGeomWindow::darea_button_signal2( GtkWidget *widget, 
					      GdkEventButton *event,
					      gpointer object )
{
    GTKGeomWindow *plotter = (GTKGeomWindow *)object;
    plotter->darea_button2( event );
    return( FALSE );
}


GTKGeomWindow::~GTKGeomWindow()
{
}


void GTKGeomWindow::update_view()
{
}


void GTKGeomWindow::zoom_fit( void )
{
    //std::cout << "Zoom fit\n";
    double min = -std::numeric_limits<double>::infinity();
    double max = std::numeric_limits<double>::infinity();
    _frame.set_ranges( PLOT_AXIS_X1, min, max );
    _frame.set_ranges( PLOT_AXIS_Y1, min, max );
    _frame.ruler_autorange_enable( PLOT_AXIS_X1, false, false );
    _frame.ruler_autorange_enable( PLOT_AXIS_Y1, false, false );

    draw_and_expose();
}


std::string GTKGeomWindow::track_text( double x, double y )
{
    double xc[3];
    int    i[3];
    std::stringstream ss;

    xc[_geomplot.vb(0)] = x;
    xc[_geomplot.vb(1)] = y;
    xc[_geomplot.vb(2)] = _geomplot.get_level_si();
    
    i[_geomplot.vb(0)] = (int)floor( (x-_geom->origo(_geomplot.vb(0)))/_geom->h() );
    i[_geomplot.vb(1)] = (int)floor( (y-_geom->origo(_geomplot.vb(1)))/_geom->h() );
    i[_geomplot.vb(2)] = _geomplot.get_level();

    Vec3D loc(xc[0],xc[1],xc[2]);

    ss << "x = " << xc[0] << " m\n"
       << "y = " << xc[1] << " m\n"
       << "z = " << xc[2] << " m\n";
    ss << "i = " << i[0] << "\n"
       << "j = " << i[1] << "\n"
       << "k = " << i[2] << "\n";
    ss << "solid = " << _geom->inside( loc ) << "\n";
    if( _epot ) {
	ss << "epot = " << (*_epot)( loc ) << "\n";
	EpotEfield efield( *_geom, *_epot );
	ss << "efield = " << efield( loc ) << "\n";
    }
    if( _bfield )
	ss << "bfield = " << (*_bfield)( loc ) << "\n";
    if( _scharge )
	ss << "scharge = " << (*_scharge)( loc ) << "\n";

    return( ss.str() );
}

struct PreferencesData {
    GtkWidget *manual_eqlines_entry;
    GtkWidget *automatic_eqlines_spin;
    GtkWidget *particlediv_spin;
    GtkWidget *scharge_field_check;
    GtkWidget *qmdiscretation_check;
    GtkWidget *meshen_check;
};


void *GTKGeomWindow::build_preferences( GtkWidget *notebook )
{
    PreferencesData *pdata = new PreferencesData;;

    GtkWidget *vbox = gtk_vbox_new( FALSE, 0 );

    // Manual eqlines
    GtkWidget *hbox = gtk_hbox_new( TRUE, 30 );
    GtkWidget *label = gtk_label_new( "Manual eqlines" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    pdata->manual_eqlines_entry = gtk_entry_new();
    std::vector<double> eqlines_manual = _geomplot.get_eqlines_manual();
    std::string s;
    for( size_t a = 0; a < eqlines_manual.size(); a++ ) {
	s += to_string( eqlines_manual[a] );
	if( a != eqlines_manual.size()-1 )
	    s += ", ";
    }
    gtk_entry_set_text( GTK_ENTRY(pdata->manual_eqlines_entry), s.c_str() );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), pdata->manual_eqlines_entry, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // Automatic eqlines
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Automatic eqlines" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    size_t eqlines_auto = _geomplot.get_eqlines_auto();
    GtkObject *automatic_eqlines_adj = gtk_adjustment_new( eqlines_auto, 0, 1000, 1, 10, 10 );
    pdata->automatic_eqlines_spin = gtk_spin_button_new( GTK_ADJUSTMENT(automatic_eqlines_adj), 1, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), pdata->automatic_eqlines_spin, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // Particle division
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Trajectory division" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    size_t particle_div = _geomplot.get_particle_div();
    GtkObject *particlediv_adj = gtk_adjustment_new( particle_div, 0, 10000, 1, 10, 10 );
    pdata->particlediv_spin = gtk_spin_button_new( GTK_ADJUSTMENT(particlediv_adj), 1, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), pdata->particlediv_spin, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // Space charge field
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Space charge field" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    pdata->scharge_field_check = gtk_check_button_new_with_label( "on/off" );
    bool scharge_field = _geomplot.get_scharge_field();
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(pdata->scharge_field_check), scharge_field );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), pdata->scharge_field_check, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // QM discretation
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Q/M discretation" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    pdata->qmdiscretation_check = gtk_check_button_new_with_label( "on/off" );
    bool qm_discretation = _geomplot.get_qm_discretation();
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(pdata->qmdiscretation_check), qm_discretation );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), pdata->qmdiscretation_check, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // Mesh
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Mesh" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    pdata->meshen_check = gtk_check_button_new_with_label( "on/off" );
    bool mesh = _geomplot.get_mesh();
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(pdata->meshen_check), mesh );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), pdata->meshen_check, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // Add notebook page
    label = gtk_label_new( "Geometry" );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox, label );

    return( (void *)pdata );
}


void GTKGeomWindow::read_preferences( GtkWidget *notebook, void *_pdata )
{
    PreferencesData *pdata = (PreferencesData *)_pdata;

    // Eqlines
    size_t eqlines_auto = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(pdata->automatic_eqlines_spin) );
    _geomplot.set_eqlines_auto( eqlines_auto );
    std::vector<double> eqlines_manual;
    char *str = (char *)gtk_entry_get_text( GTK_ENTRY(pdata->manual_eqlines_entry) );
    while( *str != '\0' ) {
	if( !isdigit( *str ) && *str != '-' && *str != '.' && *str != '+' ) {
	    str++;
	    continue;
	}
	double val = strtod( str, &str );
	std::cout << "Manual eqline: " << val << "\n";
	eqlines_manual.push_back( val );
    }
    _geomplot.set_eqlines_manual( eqlines_manual );

    // Particle division
    size_t particle_div = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(pdata->particlediv_spin) );
    _geomplot.set_particle_div( particle_div );

    // Space charge field
    bool scharge_field = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( pdata->scharge_field_check) );
    _geomplot.set_scharge_field( scharge_field );

    // QM discretation
    bool qm_discretation = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( pdata->qmdiscretation_check) );
    _geomplot.set_qm_discretation( qm_discretation );

    // Mesh
    bool mesh = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( pdata->meshen_check) );
    _geomplot.set_mesh( mesh );
}


void GTKGeomWindow::combobox( GtkComboBox *combobox )
{
    std::cout << "Combobox\n";

    int level;
    int levelmax;
    view_e view;

    switch( gtk_combo_box_get_active( GTK_COMBO_BOX(_combobox) ) ) {
    case 0:
        view = VIEW_XY;
	levelmax = _geom->size(2);
	level = levelmax/2;
        break;
    case 1:
        view = VIEW_XZ;
	levelmax = _geom->size(1);
	level = levelmax/2;
        break;
    case 2:
        view = VIEW_YZ;
	levelmax = _geom->size(0);
	level = levelmax/2;
        break;
    default:
	throw( Error( ERROR_LOCATION, "illegal combo box tag" ) );
	break;
    }

    // Set ranges to defaults on view change
    double min = -std::numeric_limits<double>::infinity();
    double max = std::numeric_limits<double>::infinity();
    _frame.set_ranges( PLOT_AXIS_X1, min, max );
    _frame.set_ranges( PLOT_AXIS_Y1, min, max );
    _frame.ruler_autorange_enable( PLOT_AXIS_X1, false, false );
    _frame.ruler_autorange_enable( PLOT_AXIS_Y1, false, false );

    // Update spinbutton
    /*
    g_signal_handlers_block_by_func( G_OBJECT(_spinbutton),
				     (void *)spinbutton_signal,
				     (gpointer)this );
    */
    g_signal_handlers_disconnect_by_func( G_OBJECT(_spinbutton),
					  (void *)spinbutton_signal,
					  (gpointer)this );
    gtk_spin_button_set_range( GTK_SPIN_BUTTON(_spinbutton), 0, levelmax-1 );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(_spinbutton), level );
    /*
    g_signal_handlers_unblock_by_func( G_OBJECT(_spinbutton),
				       (void *)spinbutton_signal,
				       (gpointer)this );
    */
    //g_signal_stop_emission_by_name( G_OBJECT(_spinbutton), "value-changed" );
    g_signal_connect( G_OBJECT(_spinbutton), "value-changed",
                      G_CALLBACK(spinbutton_signal),
                      (gpointer)this );


    // Update view
    _geomplot.set_view( view, level );
    update_view();
    draw_and_expose();
}


void GTKGeomWindow::spinbutton( GtkSpinButton *spinbutton )
{
    std::cout << "Spinbutton\n";
    //g_signal_stop_emission_by_name( G_OBJECT(_spinbutton), "value-changed" );

    int level = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(_spinbutton) );

    std::stringstream ss;
    switch( _geomplot.get_view() ) {
    case VIEW_XY:
	ss << "z = "<< _geom->origo(2)+level*_geom->h() << " m";
	break;
    case VIEW_XZ:
	ss << "y = "<< _geom->origo(1)+level*_geom->h() << " m";
	break;
    case VIEW_YZ:
	ss << "x = "<< _geom->origo(0)+level*_geom->h() << " m";
	break;
    }
    gtk_statusbar_pop( GTK_STATUSBAR(_statusbar), 0 );
    gtk_statusbar_push( GTK_STATUSBAR(_statusbar), 0, ss.str().c_str() );

    // Update view
    _geomplot.set_view( _geomplot.get_view(), level );
    update_view();
    draw_and_expose();
}


void GTKGeomWindow::menuitem_tool_change( GtkToolButton *button )
{
    int tool;
    const char *label = gtk_tool_button_get_label( button );
    if( !strcmp( label, "Particle diagnostics" ) ) 
	tool = TOOL_PARTICLE_DIAG;
    else if( !strcmp( label, "Field diagnostics" ) ) 
	tool = TOOL_FIELD_DIAG;
    else {
	tool = TOOL_UNKNOWN;
	return;
    }

    if( !gtk_toggle_tool_button_get_active( GTK_TOGGLE_TOOL_BUTTON(button) ) ) {
	// Disable tool
	_tool = TOOL_UNKNOWN;
    } else {
	// Enable tool
	_tool = tool;
    }
}


void GTKGeomWindow::combobox_signal( GtkComboBox *combobox,
				     gpointer object )
{
    GTKGeomWindow *plotter = (GTKGeomWindow *)object;
    plotter->combobox( combobox );
}


void GTKGeomWindow::spinbutton_signal( GtkSpinButton *spinbutton,
				       gpointer object )
{
    GTKGeomWindow *plotter = (GTKGeomWindow *)object;
    plotter->spinbutton( spinbutton );
}



void GTKGeomWindow::menuitem_tool_change_signal( GtkToolButton *button,
						 gpointer object )
{
    GTKGeomWindow *plotter = (GTKGeomWindow *)object;
    plotter->menuitem_tool_change( button );
}

