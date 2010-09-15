/*! \file gtkparticlediagdialog.cpp
 *  \brief Source code for gtkparticlediagdialog.cpp
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

#include "gtkparticlediagdialog.hpp"
#include "gtkparticlediagwindow.hpp"
#include "types.hpp"


GTKParticleDiagDialog::GTKParticleDiagDialog( GtkWidget *window, GTKPlotter *plotter, int plane, double val )
    : _window(window), _plotter(plotter), _plane(plane), _val(val)
{
    _geom = _plotter->get_geometry();
}


GTKParticleDiagDialog::~GTKParticleDiagDialog()
{

}


void GTKParticleDiagDialog::profile_toggled2( GtkToggleButton *togglebutton )
{
    if( gtk_toggle_button_get_active( togglebutton ) ) {
	gtk_widget_set_sensitive( _radio_plot_scatter, FALSE );
	gtk_widget_set_sensitive( _radio_plot_colormap, FALSE );
    } else {
	gtk_widget_set_sensitive( _radio_plot_scatter, TRUE );
	gtk_widget_set_sensitive( _radio_plot_colormap, TRUE );
    }
}

void GTKParticleDiagDialog::profile_toggled( GtkToggleButton *togglebutton,
					     gpointer         user_data )
{
    GTKParticleDiagDialog *dialog = (GTKParticleDiagDialog *)user_data;
    dialog->profile_toggled2( togglebutton );
}


void GTKParticleDiagDialog::plane_activated( void )
{
    // Plane changed
    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_plane_x) ) ) {
	gtk_widget_set_sensitive( _radio_emit_xx, FALSE );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_emit_xx) ) )
	    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_yy), TRUE );
	gtk_widget_set_sensitive( _radio_emit_yy, TRUE );
	if( _geom->geom_mode() == MODE_3D ) {
	    gtk_widget_set_sensitive( _radio_emit_zz, TRUE );
	    gtk_widget_set_sensitive( _radio_prof_yz, TRUE );
	    gtk_widget_set_sensitive( _radio_prof_xz, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xz) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_yy), TRUE );
	    gtk_widget_set_sensitive( _radio_prof_xy, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xy) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_yy), TRUE );
	} else {
	    gtk_widget_set_sensitive( _radio_prof_yz, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_yz) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_yy), TRUE );
	    gtk_widget_set_sensitive( _radio_prof_xz, TRUE );
	}


    } else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_plane_y) ) ) {
	gtk_widget_set_sensitive( _radio_emit_xx, TRUE );
	gtk_widget_set_sensitive( _radio_emit_yy, FALSE );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_emit_yy) ) )
	    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	if( _geom->geom_mode() == MODE_3D ) {
	    gtk_widget_set_sensitive( _radio_emit_zz, TRUE );
	    gtk_widget_set_sensitive( _radio_prof_yz, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_yz) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	    gtk_widget_set_sensitive( _radio_prof_xz, TRUE );
	    gtk_widget_set_sensitive( _radio_prof_xy, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xy) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );

	} else {
	    gtk_widget_set_sensitive( _radio_prof_yz, TRUE );
	    gtk_widget_set_sensitive( _radio_prof_xz, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xz) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	}

    } else {
	gtk_widget_set_sensitive( _radio_emit_xx, TRUE );
	gtk_widget_set_sensitive( _radio_emit_yy, TRUE );
	if( _geom->geom_mode() == MODE_3D ) {
	    gtk_widget_set_sensitive( _radio_emit_zz, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_emit_zz) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	    gtk_widget_set_sensitive( _radio_prof_yz, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_yz) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	    gtk_widget_set_sensitive( _radio_prof_xz, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xz) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	    gtk_widget_set_sensitive( _radio_prof_xy, TRUE );
	}
    }
}


void GTKParticleDiagDialog::plane_toggled( GtkToggleButton *togglebutton,
					   gpointer         user_data )
{
    GTKParticleDiagDialog *dialog = (GTKParticleDiagDialog *)user_data;
    if( gtk_toggle_button_get_active( togglebutton ) )
	dialog->plane_activated();
}


void GTKParticleDiagDialog::run( void )
{
    GtkWidget *dialog = gtk_dialog_new_with_buttons( "Make particle diagnostics",
						     GTK_WINDOW(_window),
						     (GtkDialogFlags)(GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT), 
						     GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
						     GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
						     NULL );
    GtkWidget *vbox = GTK_DIALOG(dialog)->vbox;

    // ****************************************************************************

    GtkWidget *hbox = gtk_hbox_new( TRUE, 30 );
    GtkWidget *vbox2 = gtk_vbox_new( FALSE, 0 );

    // Label for plane
    GtkWidget *label = gtk_label_new( "" );
    gtk_label_set_markup( GTK_LABEL(label), "<span weight=\"bold\">Plane</span>" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0 );
    gtk_box_pack_start( GTK_BOX(vbox2), label, FALSE, TRUE, 0 );

    // Radio buttons for selecting plane
    GtkWidget *alignment = gtk_alignment_new( 0, 0, 0, 0 );
    gtk_alignment_set_padding( GTK_ALIGNMENT(alignment), 0, 0, 15, 0 );
    GtkWidget *vbox3 = gtk_vbox_new( FALSE, 0 );
    _radio_plane_x = gtk_radio_button_new_with_label_from_widget( NULL,
								  "X-axis" );
    if( _plane == 0 ) gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_plane_x), TRUE );
    else  gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_plane_x), FALSE );
    g_signal_connect( _radio_plane_x, "toggled",
		      G_CALLBACK(plane_toggled), (gpointer)this );
    gtk_box_pack_start( GTK_BOX(vbox3), _radio_plane_x, FALSE, TRUE, 0 );

    if( _geom->geom_mode() == MODE_2D || _geom->geom_mode() == MODE_3D )
	_radio_plane_y = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_plane_x),
								      "Y-axis");
    else
	_radio_plane_y = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_plane_x),
								      "R-axis");
    if( _plane == 1 ) gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_plane_y), TRUE );
    else  gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_plane_y), FALSE );
    g_signal_connect( _radio_plane_y, "toggled",
		      G_CALLBACK(plane_toggled), (gpointer)this );
    gtk_box_pack_start( GTK_BOX(vbox3), _radio_plane_y, FALSE, TRUE, 0 );

    if( _geom->geom_mode() == MODE_3D ) {
	_radio_plane_z = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_plane_x),
								      "Z-axis");
	if( _plane == 2 ) gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_plane_z), TRUE );
	else  gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_plane_z), FALSE );
	g_signal_connect( _radio_plane_z, "toggled",
			  G_CALLBACK(plane_toggled), (gpointer)this );
	gtk_box_pack_start( GTK_BOX(vbox3), _radio_plane_z, FALSE, TRUE, 0 );
    }

    gtk_container_add( GTK_CONTAINER(alignment), vbox3 );
    gtk_box_pack_start( GTK_BOX(vbox2), alignment, FALSE, TRUE, 0 );

    gtk_box_pack_start( GTK_BOX(hbox), vbox2, FALSE, TRUE, 0 );

    // ****************************************************************************

    vbox2 = gtk_vbox_new( FALSE, 0 );

    // Label for level
    label = gtk_label_new( "" );
    gtk_label_set_markup( GTK_LABEL(label), "<span weight=\"bold\">Level</span>" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0 );
    gtk_box_pack_start( GTK_BOX(vbox2), label, FALSE, TRUE, 0 );

    // Level
    alignment = gtk_alignment_new( 0, 0, 0, 0 );
    gtk_alignment_set_padding( GTK_ALIGNMENT(alignment), 0, 0, 15, 0 );
    GtkWidget *entry_level = gtk_entry_new_with_max_length( 30 );
    char buf[128];
    snprintf( buf, 128, "%lf", _val );
    gtk_entry_set_text( GTK_ENTRY(entry_level), buf );
    gtk_container_add( GTK_CONTAINER(alignment), entry_level );
    gtk_box_pack_start( GTK_BOX(vbox2), alignment, FALSE, TRUE, 0 );

    gtk_box_pack_start( GTK_BOX(hbox), vbox2, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // ****************************************************************************

    // Separator
    GtkWidget *separator = gtk_hseparator_new();
    gtk_box_pack_start( GTK_BOX(vbox), separator, FALSE, TRUE, 0 );

    // Label for plot
    label = gtk_label_new( "" );
    gtk_label_set_markup( GTK_LABEL(label), "<span weight=\"bold\">Plot</span>" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), label, FALSE, TRUE, 0 );

    // Columns for radio buttons
    hbox = gtk_hbox_new( TRUE, 30 );

    // Radio buttons for plot type
    alignment = gtk_alignment_new( 0, 0, 0, 0 );
    gtk_alignment_set_padding( GTK_ALIGNMENT(alignment), 0, 0, 15, 0 );
    vbox2 = gtk_vbox_new( FALSE, 0 );
    _radio_emit_xx = gtk_radio_button_new_with_label_from_widget( NULL,
								  "Emittance x-x\'" );
    if( _plane != 0 )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
    gtk_box_pack_start( GTK_BOX(vbox2), _radio_emit_xx, FALSE, TRUE, 0 );
    if( _geom->geom_mode() == MODE_CYL )
	_radio_emit_yy = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								      "Emittance r-r\'" );
    else
	_radio_emit_yy = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								      "Emittance y-y\'" );
    if( _plane == 0 )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_yy), TRUE );
    gtk_box_pack_start( GTK_BOX(vbox2), _radio_emit_yy, FALSE, TRUE, 0 );
    if( _geom->geom_mode() == MODE_3D ) {
	_radio_emit_zz = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								      "Emittance z-z\'" );
	gtk_box_pack_start( GTK_BOX(vbox2), _radio_emit_zz, FALSE, TRUE, 0 );
    }
    if( _geom->geom_mode() == MODE_3D ) {
	_radio_prof_yz = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								      "Profile y-z" );
	gtk_box_pack_start( GTK_BOX(vbox2), _radio_prof_yz, FALSE, TRUE, 0 );
	_radio_prof_xz = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								      "Profile x-z" );
	gtk_box_pack_start( GTK_BOX(vbox2), _radio_prof_xz, FALSE, TRUE, 0 );
	_radio_prof_xy = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								      "Profile x-y" );
	gtk_box_pack_start( GTK_BOX(vbox2), _radio_prof_xy, FALSE, TRUE, 0 );
    } else {
	_radio_prof_yz = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								      "Profile x" );
	g_signal_connect( _radio_prof_yz, "toggled",
			  G_CALLBACK(profile_toggled), (gpointer)this );
	gtk_box_pack_start( GTK_BOX(vbox2), _radio_prof_yz, FALSE, TRUE, 0 );
	if( _geom->geom_mode() == MODE_2D )
	    _radio_prof_xz = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
									  "Profile y" );
	else
	    _radio_prof_xz = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
									  "Profile r" );
	g_signal_connect( _radio_prof_xz, "toggled",
			  G_CALLBACK(profile_toggled), (gpointer)this );
	gtk_box_pack_start( GTK_BOX(vbox2), _radio_prof_xz, FALSE, TRUE, 0 );
    }
    gtk_container_add( GTK_CONTAINER(alignment), vbox2 );
    gtk_box_pack_start( GTK_BOX(hbox), alignment, FALSE, TRUE, 0 );

    // Radio buttons for plot style
    alignment = gtk_alignment_new( 0, 0, 0, 0 );
    gtk_alignment_set_padding( GTK_ALIGNMENT(alignment), 0, 0, 15, 0 );
    vbox2 = gtk_vbox_new( FALSE, 0 );
    _radio_plot_scatter = gtk_radio_button_new_with_label_from_widget( NULL,
								       "Scatter" );
    gtk_box_pack_start( GTK_BOX(vbox2), _radio_plot_scatter, FALSE, TRUE, 0 );
    _radio_plot_colormap = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_plot_scatter),
									"Colormap" );
    gtk_box_pack_start( GTK_BOX(vbox2), _radio_plot_colormap, FALSE, TRUE, 0 );
    gtk_container_add( GTK_CONTAINER(alignment), vbox2 );
    gtk_box_pack_start( GTK_BOX(hbox), alignment, FALSE, TRUE, 0 );

    // End
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // ****************************************************************************

    plane_activated();

    gtk_widget_show_all( dialog );
    particle_diag_plot_type_e type = PARTICLE_DIAG_PLOT_SCATTER;
    if( gtk_dialog_run( GTK_DIALOG(dialog) ) == GTK_RESPONSE_ACCEPT ) {

	// Read in options and start particle diagnostics plot

	// Read plane
	coordinate_axis_e axis;
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_plane_x) ) )
	    axis = AXIS_X;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_plane_y) ) )
	    axis = AXIS_Y;
	else if( _geom->geom_mode() == MODE_CYL && 
		 gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_plane_y) ) )
	    axis = AXIS_R;
	else
	    axis = AXIS_Z;
	//std::cout << "axis = " << axis << "\n";

	// Read level
	double level;
	const char *entry_text = gtk_entry_get_text( GTK_ENTRY(entry_level) );
	level = atof( entry_text );
	//std::cout << "level = " << level << "\n";

	// Read plot type
	trajectory_diagnostic_e diagx, diagy = DIAG_XP;
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_emit_xx) ) ) {
	    diagx = DIAG_X;
	    diagy = DIAG_XP;
	} else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_emit_yy) ) ) {
	    if( _geom->geom_mode() == MODE_CYL ) {
		diagx = DIAG_R;
		diagy = DIAG_RP;
	    } else {
		diagx = DIAG_Y;
		diagy = DIAG_YP;
	    }
	} else if( _geom->geom_mode() == MODE_3D && 
		   gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_emit_zz) ) ) {
	    diagx = DIAG_Z;
	    diagy = DIAG_ZP;
	} else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_yz) ) ) {
	    if( _geom->geom_mode() == MODE_3D ) {
		diagx = DIAG_Y;
		diagy = DIAG_Z;
	    } else {
		diagx = DIAG_X;
		type = PARTICLE_DIAG_PLOT_HISTO1D;
	    }
	} else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xz) ) ) {
	    if( _geom->geom_mode() == MODE_3D ) {
		diagx = DIAG_X;
		diagy = DIAG_Z;
	    } else if( _geom->geom_mode() == MODE_2D ) {
		diagx = DIAG_Y;
		type = PARTICLE_DIAG_PLOT_HISTO1D;
	    } else {
		diagx = DIAG_R;
		type = PARTICLE_DIAG_PLOT_HISTO1D;
	    }
	} else {
	    diagx = DIAG_X;
	    diagy = DIAG_Y;
	}

	// Read plot style
	if( type != PARTICLE_DIAG_PLOT_HISTO1D ) {
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_plot_scatter) ) )
		type = PARTICLE_DIAG_PLOT_SCATTER;
	    else
		type = PARTICLE_DIAG_PLOT_HISTO2D;
	}

	_plotter->new_particle_plot_window( axis, level, type, diagx, diagy );
    }

    gtk_widget_destroy( dialog );
}















