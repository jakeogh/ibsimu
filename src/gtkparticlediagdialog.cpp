/*! \file gtkparticlediagdialog.cpp
 *  \brief Source code for gtkparticlediagdialog.cpp
 */

/* Copyright (c) 2005-2012 Taneli Kalvas. All rights reserved.
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

#include "gtkparticlediagdialog.hpp"
#include "gtkparticlediagwindow.hpp"
#include "types.hpp"


GTKParticleDiagDialog::GTKParticleDiagDialog( GtkWidget *window, GTKPlotter &plotter, 
					      int plane, double val )
    : _window(window), _plotter(&plotter), _plane(plane), _val(val),
    _radio_plane_x(NULL), _radio_plane_y(NULL), _radio_plane_z(NULL),
    _radio_emit_xx(NULL), _radio_emit_yy(NULL), _radio_emit_zz(NULL),
    _radio_prof_yz(NULL), _radio_prof_xz(NULL), _radio_prof_xy(NULL),
    _radio_plot_scatter(NULL), _radio_plot_colormap(NULL),
    _radio_prof_x(NULL), _radio_prof_y(NULL), _radio_prof_z(NULL),
    _radio_prof_xp(NULL), _radio_prof_yp(NULL), _radio_prof_zp(NULL),
    _radio_energy(NULL), _radio_qm(NULL), _radio_charge(NULL), _radio_mass(NULL)
{
    _geom = _plotter->get_geometry();
}


GTKParticleDiagDialog::~GTKParticleDiagDialog()
{

}


void GTKParticleDiagDialog::plot1d_toggled2( GtkToggleButton *togglebutton )
{
    if( gtk_toggle_button_get_active( togglebutton ) ) {
	gtk_widget_set_sensitive( _radio_plot_scatter, FALSE );
	gtk_widget_set_sensitive( _radio_plot_colormap, FALSE );
    } else {
	gtk_widget_set_sensitive( _radio_plot_scatter, TRUE );
	gtk_widget_set_sensitive( _radio_plot_colormap, TRUE );
    }
}


void GTKParticleDiagDialog::plot1d_toggled( GtkToggleButton *togglebutton,
					     gpointer         user_data )
{
    GTKParticleDiagDialog *dialog = (GTKParticleDiagDialog *)user_data;
    dialog->plot1d_toggled2( togglebutton );
}


void GTKParticleDiagDialog::conversion_toggled2( GtkToggleButton *togglebutton )
{
    if( gtk_toggle_button_get_active( togglebutton ) ) {
	gtk_widget_set_sensitive( _radio_plot_scatter, FALSE );
	gtk_widget_set_sensitive( _radio_plot_colormap, TRUE );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_plot_colormap), TRUE );
    } else {
	gtk_widget_set_sensitive( _radio_plot_scatter, TRUE );
	gtk_widget_set_sensitive( _radio_plot_colormap, TRUE );
    }
}


void GTKParticleDiagDialog::conversion_toggled( GtkToggleButton *togglebutton,
						gpointer         user_data )
{
    GTKParticleDiagDialog *dialog = (GTKParticleDiagDialog *)user_data;
    dialog->conversion_toggled2( togglebutton );
}


void GTKParticleDiagDialog::plane_activated( void )
{
    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_plane_x) ) ) {

	gtk_widget_set_sensitive( _radio_prof_x, FALSE );
	gtk_widget_set_sensitive( _radio_prof_xp, FALSE );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_x) ) ||
	    gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xp) ) )
	    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_yy), TRUE );

	gtk_widget_set_sensitive( _radio_emit_xx, FALSE );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_emit_xx) ) )
	    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_yy), TRUE );

	gtk_widget_set_sensitive( _radio_prof_y, TRUE );
	gtk_widget_set_sensitive( _radio_prof_yp, TRUE );
	gtk_widget_set_sensitive( _radio_emit_yy, TRUE );

	if( _geom->geom_mode() == MODE_CYL ) {
	    // EmittanceConv
	    gtk_widget_set_sensitive( _radio_emit_zz, TRUE );
	}
	if( _geom->geom_mode() == MODE_3D ) {
	    gtk_widget_set_sensitive( _radio_prof_z, TRUE );
	    gtk_widget_set_sensitive( _radio_prof_zp, TRUE );
	    gtk_widget_set_sensitive( _radio_emit_zz, TRUE );
	    gtk_widget_set_sensitive( _radio_prof_yz, TRUE );
	    gtk_widget_set_sensitive( _radio_prof_xz, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xz) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_yy), TRUE );
	    gtk_widget_set_sensitive( _radio_prof_xy, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xy) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_yy), TRUE );
	} else {
	    gtk_widget_set_sensitive( _radio_prof_x, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_x) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_yy), TRUE );
	    gtk_widget_set_sensitive( _radio_prof_y, TRUE );
	}


    } else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_plane_y) ) ) {

	gtk_widget_set_sensitive( _radio_prof_x, TRUE );
	gtk_widget_set_sensitive( _radio_prof_xp, TRUE );
	gtk_widget_set_sensitive( _radio_emit_xx, TRUE );

	gtk_widget_set_sensitive( _radio_prof_y, FALSE );
	gtk_widget_set_sensitive( _radio_prof_yp, FALSE );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_y) ) ||
	    gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_yp) ) )
	    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );

	gtk_widget_set_sensitive( _radio_emit_yy, FALSE );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_emit_yy) ) )
	    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );

	if( _geom->geom_mode() == MODE_CYL ) {
	    // EmittanceConv
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_emit_zz) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	    gtk_widget_set_sensitive( _radio_emit_zz, FALSE );
	}
	if( _geom->geom_mode() == MODE_3D ) {
	    gtk_widget_set_sensitive( _radio_prof_z, TRUE );
	    gtk_widget_set_sensitive( _radio_prof_zp, TRUE );
	    gtk_widget_set_sensitive( _radio_emit_zz, TRUE );
	    gtk_widget_set_sensitive( _radio_prof_yz, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_yz) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	    gtk_widget_set_sensitive( _radio_prof_xz, TRUE );
	    gtk_widget_set_sensitive( _radio_prof_xy, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xy) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );

	} else {
	    gtk_widget_set_sensitive( _radio_prof_x, TRUE );
	    gtk_widget_set_sensitive( _radio_prof_y, FALSE );
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_y) ) )
		gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	}

    } else /* gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_plane_z) ) */ {

	gtk_widget_set_sensitive( _radio_prof_x, TRUE );
	gtk_widget_set_sensitive( _radio_prof_xp, TRUE );
	gtk_widget_set_sensitive( _radio_emit_xx, TRUE );

	gtk_widget_set_sensitive( _radio_prof_y, TRUE );
	gtk_widget_set_sensitive( _radio_prof_yp, TRUE );
	gtk_widget_set_sensitive( _radio_emit_yy, TRUE );

	gtk_widget_set_sensitive( _radio_emit_zz, FALSE );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_emit_zz) ) )
	    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	
	gtk_widget_set_sensitive( _radio_prof_yz, FALSE );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_yz) ) )
	    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	
	gtk_widget_set_sensitive( _radio_prof_xz, FALSE );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xz) ) )
	    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	
	gtk_widget_set_sensitive( _radio_prof_z, FALSE );
	gtk_widget_set_sensitive( _radio_prof_zp, FALSE );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_z) ) ||
	    gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_zp) ) )
	    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_xx), TRUE );
	
	gtk_widget_set_sensitive( _radio_prof_xy, TRUE );
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
						     (GtkDialogFlags)(GTK_DIALOG_MODAL | 
								      GTK_DIALOG_DESTROY_WITH_PARENT), 
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
    snprintf( buf, 128, "%g", _val );
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
    gtk_label_set_markup( GTK_LABEL(label), "<span weight=\"bold\">Plot 2D</span>" );
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
    if( _geom->geom_mode() == MODE_CYL ) {
	_radio_emit_yy = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								      "Emittance r-r\'" );
	if( _plane == 0 )
	    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_yy), TRUE );
	gtk_box_pack_start( GTK_BOX(vbox2), _radio_emit_yy, FALSE, TRUE, 0 );
	_radio_emit_zz = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								      "Emittance z-z\' (conv)" );
	g_signal_connect( _radio_emit_zz, "toggled",
			  G_CALLBACK(conversion_toggled), (gpointer)this );
	gtk_box_pack_start( GTK_BOX(vbox2), _radio_emit_zz, FALSE, TRUE, 0 );
    } else {
	_radio_emit_yy = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								      "Emittance y-y\'" );
	if( _plane == 0 )
	    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_radio_emit_yy), TRUE );
	gtk_box_pack_start( GTK_BOX(vbox2), _radio_emit_yy, FALSE, TRUE, 0 );
    }
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
    }
    gtk_container_add( GTK_CONTAINER(alignment), vbox2 );
    gtk_box_pack_start( GTK_BOX(hbox), alignment, FALSE, TRUE, 0 );

    // Radio buttons for 2d plot style
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

    // Separator
    separator = gtk_hseparator_new();
    gtk_box_pack_start( GTK_BOX(vbox), separator, FALSE, TRUE, 0 );

    // Label for plot
    label = gtk_label_new( "" );
    gtk_label_set_markup( GTK_LABEL(label), "<span weight=\"bold\">Plot 1D</span>" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), label, FALSE, TRUE, 0 );

    // Columns for radio buttons
    hbox = gtk_hbox_new( TRUE, 30 );

    // First column
    alignment = gtk_alignment_new( 0, 0, 0, 0 );
    gtk_alignment_set_padding( GTK_ALIGNMENT(alignment), 0, 0, 15, 0 );
    vbox2 = gtk_vbox_new( FALSE, 0 );

    // x
    _radio_prof_x = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								 "Profile x" );
    g_signal_connect( _radio_prof_x, "toggled",
		      G_CALLBACK(plot1d_toggled), (gpointer)this );
    gtk_box_pack_start( GTK_BOX(vbox2), _radio_prof_x, FALSE, TRUE, 0 );

    // y
    if( _geom->geom_mode() == MODE_CYL )
	_radio_prof_y = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								     "Profile r" );
    else
	_radio_prof_y = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								     "Profile y" );
    g_signal_connect( _radio_prof_y, "toggled",
		      G_CALLBACK(plot1d_toggled), (gpointer)this );
    gtk_box_pack_start( GTK_BOX(vbox2), _radio_prof_y, FALSE, TRUE, 0 );

    // z
    if( _geom->geom_mode() == MODE_3D ) {
	_radio_prof_z = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								 "Profile z" );
	g_signal_connect( _radio_prof_z, "toggled",
			  G_CALLBACK(plot1d_toggled), (gpointer)this );
	gtk_box_pack_start( GTK_BOX(vbox2), _radio_prof_z, FALSE, TRUE, 0 );
    }

    // energy
    _radio_energy = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								 "Energy" );
    //g_signal_connect( _radio_energy, "toggled",
    //G_CALLBACK(plot1d_toggled), (gpointer)this );
    gtk_box_pack_start( GTK_BOX(vbox2), _radio_energy, FALSE, TRUE, 0 );

    gtk_container_add( GTK_CONTAINER(alignment), vbox2 );
    gtk_box_pack_start( GTK_BOX(hbox), alignment, FALSE, TRUE, 0 );

    // charge
    _radio_charge = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								 "Charge" );
    //g_signal_connect( _radio_charge, "toggled",
    //G_CALLBACK(plot1d_toggled), (gpointer)this );
    gtk_box_pack_start( GTK_BOX(vbox2), _radio_charge, FALSE, TRUE, 0 );


    // Second column
    alignment = gtk_alignment_new( 0, 0, 0, 0 );
    gtk_alignment_set_padding( GTK_ALIGNMENT(alignment), 0, 0, 15, 0 );
    vbox2 = gtk_vbox_new( FALSE, 0 );

    // x'
    _radio_prof_xp = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								 "Profile x'" );
    g_signal_connect( _radio_prof_xp, "toggled",
		      G_CALLBACK(plot1d_toggled), (gpointer)this );
    gtk_box_pack_start( GTK_BOX(vbox2), _radio_prof_xp, FALSE, TRUE, 0 );

    // y'
    if( _geom->geom_mode() == MODE_CYL )
	_radio_prof_yp = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								     "Profile r'" );
    else
	_radio_prof_yp = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								     "Profile y'" );
    g_signal_connect( _radio_prof_yp, "toggled",
		      G_CALLBACK(plot1d_toggled), (gpointer)this );
    gtk_box_pack_start( GTK_BOX(vbox2), _radio_prof_yp, FALSE, TRUE, 0 );

    // z'
    if( _geom->geom_mode() == MODE_3D ) {
	_radio_prof_zp = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
								 "Profile z'" );
	g_signal_connect( _radio_prof_zp, "toggled",
			  G_CALLBACK(plot1d_toggled), (gpointer)this );
	gtk_box_pack_start( GTK_BOX(vbox2), _radio_prof_zp, FALSE, TRUE, 0 );
    }

    // q/m
    _radio_qm = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
							     "Q/M" );
    //g_signal_connect( _radio_qm, "toggled",
    //G_CALLBACK(plot1d_toggled), (gpointer)this );
    gtk_box_pack_start( GTK_BOX(vbox2), _radio_qm, FALSE, TRUE, 0 );

    // mass
    _radio_mass = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_radio_emit_xx),
							       "Mass" );
    //g_signal_connect( _radio_mass, "toggled",
    //G_CALLBACK(plot1d_toggled), (gpointer)this );
    gtk_box_pack_start( GTK_BOX(vbox2), _radio_mass, FALSE, TRUE, 0 );

    gtk_container_add( GTK_CONTAINER(alignment), vbox2 );
    gtk_box_pack_start( GTK_BOX(hbox), alignment, FALSE, TRUE, 0 );

    // End
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // ****************************************************************************

    // Set active/passive options
    plane_activated();

    gtk_widget_show_all( dialog );
    particle_diag_plot_type_e type = PARTICLE_DIAG_PLOT_NONE;
    if( gtk_dialog_run( GTK_DIALOG(dialog) ) == GTK_RESPONSE_ACCEPT ) {

	// Read in options and start particle diagnostics plot

	// Read plane
	coordinate_axis_e axis;
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_plane_x) ) )
	    axis = AXIS_X;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_plane_y) ) ) {
	    if( _geom->geom_mode() == MODE_CYL )
		axis = AXIS_R;
	    else
		axis = AXIS_Y;
	} else
	    axis = AXIS_Z;
	//std::cout << "axis = " << axis << "\n";

	// Read level
	double level;
	const char *entry_text = gtk_entry_get_text( GTK_ENTRY(entry_level) );
	level = atof( entry_text );
	//std::cout << "level = " << level << "\n";

	// Read plot type and set diagnostic accordingly
	trajectory_diagnostic_e diagx = DIAG_NONE, diagy = DIAG_NONE;
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
	} else if( _radio_emit_zz &&
		   gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_emit_zz) ) ) {
	    if( _geom->geom_mode() == MODE_3D ) {
		diagx = DIAG_Z;
		diagy = DIAG_ZP;
	    } else /* MODE_CYL */ {
		// EmittanceConv
		diagx = DIAG_Z;
		diagy = DIAG_ZP;
		type = PARTICLE_DIAG_PLOT_HISTO2D;
	    }
	} else if( _radio_prof_yz &&
		   gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_yz) ) ) {
	    diagx = DIAG_Y;
	    diagy = DIAG_Z;
	} else if( _radio_prof_xz &&
		   gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xz) ) ) {
	    diagx = DIAG_X;
	    diagy = DIAG_Z;
	} else if( _radio_prof_xy &&
		   gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xy) ) ) {
	    diagx = DIAG_X;
	    diagy = DIAG_Y;
	} else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_x) ) ) {
	    diagx = DIAG_X;
	    type = PARTICLE_DIAG_PLOT_HISTO1D;
	} else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_y) ) ) {
	    if( _geom->geom_mode() == MODE_CYL )
		diagx = DIAG_R;
	    else
		diagx = DIAG_Y;
	    type = PARTICLE_DIAG_PLOT_HISTO1D;
	} else if( _radio_prof_z && 
		   gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_z) ) ) {
	    diagx = DIAG_Z;
	    type = PARTICLE_DIAG_PLOT_HISTO1D;
	} else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_xp) ) ) {
	    diagx = DIAG_XP;
	    type = PARTICLE_DIAG_PLOT_HISTO1D;
	} else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_yp) ) ) {
	    if( _geom->geom_mode() == MODE_CYL )
		diagx = DIAG_RP;
	    else
		diagx = DIAG_YP;
	    type = PARTICLE_DIAG_PLOT_HISTO1D;
	} else if( _radio_prof_zp &&
		   gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_prof_zp) ) ) {
	    diagx = DIAG_ZP;
	    type = PARTICLE_DIAG_PLOT_HISTO1D;
	} else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_energy) ) ) {
	    diagx = DIAG_EK;
	    type = PARTICLE_DIAG_PLOT_HISTO1D;
	} else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_qm) ) ) {
	    diagx = DIAG_QM;
	    type = PARTICLE_DIAG_PLOT_HISTO1D;
	} else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_mass) ) ) {
	    diagx = DIAG_MASS;
	    type = PARTICLE_DIAG_PLOT_HISTO1D;
	} else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_charge) ) ) {
	    diagx = DIAG_CHARGE;
	    type = PARTICLE_DIAG_PLOT_HISTO1D;
	} else {
	    throw( Error( ERROR_LOCATION, "unknown diagnostic (internal error)" ) );
	}

	// Read 2d plot style if not set already
	if( type == PARTICLE_DIAG_PLOT_NONE ) {
	    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_radio_plot_scatter) ) )
		type = PARTICLE_DIAG_PLOT_SCATTER;
	    else
		type = PARTICLE_DIAG_PLOT_HISTO2D;
	}

	_plotter->new_particle_plot_window( axis, level, type, diagx, diagy );
    }

    gtk_widget_destroy( dialog );
}








