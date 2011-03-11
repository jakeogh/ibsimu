/*! \file gtkparticlediagwindow.cpp
 *  \brief %Particle diagnostic window.
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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

#include <sstream>
#include <limits>
#include "gtkparticlediagwindow.hpp"
#include "gtkparticlediagexportdialog.hpp"
#include "histogram.hpp"


GTKParticleDiagWindow::GTKParticleDiagWindow( GTKPlotter *plotter, const ParticleDataBase *pdb, 
					      const Geometry *geom,
					      coordinate_axis_e axis, double level, 
					      particle_diag_plot_type_e type,
					      trajectory_diagnostic_e diagx, 
					      trajectory_diagnostic_e diagy )
    : GTKWindow(plotter), _plot(&_frame, geom, pdb, axis, level, type, diagx, diagy)
{
    // Set window title
    gtk_window_set_title( GTK_WINDOW(_window), "Particle diagnostics" );

    // Add export menu item
    GtkWidget *item_export = gtk_menu_item_new_with_mnemonic( "_Export" );
    gtk_menu_shell_prepend( GTK_MENU_SHELL(_menu_file), item_export );
    g_signal_connect( G_OBJECT(item_export), "activate",
                      G_CALLBACK(menuitem_export_signal),
                      (gpointer)this );

    _plot.build_plot();
    show();
}


GTKParticleDiagWindow::~GTKParticleDiagWindow()
{
}



struct PreferencesData {
    GtkWidget *radio_plot_scatter;
    GtkWidget *radio_plot_colormap;
    GtkWidget *histo_n_spin;
    GtkWidget *histo_m_spin;
    GtkWidget *radio_int_closest;
    GtkWidget *radio_int_bilinear;
    GtkWidget *radio_int_bicubic;
    GtkWidget *dot_size_spin;
    GtkWidget *ellipse_check;
};


void *GTKParticleDiagWindow::build_preferences( GtkWidget *notebook )
{
    PreferencesData *pdata = new PreferencesData;;

    GtkWidget *vbox = gtk_vbox_new( FALSE, 0 );

    // Plot type
    GtkWidget *hbox = gtk_hbox_new( TRUE, 30 );
    GtkWidget *label = gtk_label_new( "Plot type" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    GtkWidget *vbox2 = gtk_vbox_new( FALSE, 0 );
    pdata->radio_plot_scatter = gtk_radio_button_new_with_label_from_widget( NULL,
								       "Scatter" );

    // Get current parameters
    particle_diag_plot_type_e type = _plot.get_type();
    interpolation_e interpolation =  _plot.get_colormap_interpolation();

    if( type == PARTICLE_DIAG_PLOT_SCATTER )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(pdata->radio_plot_scatter), true );
    gtk_box_pack_start( GTK_BOX(vbox2), pdata->radio_plot_scatter, FALSE, TRUE, 0 );
    pdata->radio_plot_colormap = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(pdata->radio_plot_scatter),
									       "Colormap" );
    if( type == PARTICLE_DIAG_PLOT_HISTO2D )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(pdata->radio_plot_colormap), true );
    gtk_box_pack_start( GTK_BOX(vbox2), pdata->radio_plot_colormap, FALSE, TRUE, 0 );
    if( type == PARTICLE_DIAG_PLOT_HISTO1D ) {
	gtk_widget_set_sensitive( pdata->radio_plot_scatter, false );
	gtk_widget_set_sensitive( pdata->radio_plot_colormap, false );
    }
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), vbox2, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // Histogram size
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Bin size n" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    GtkObject *histo_n_adj = gtk_adjustment_new( _plot.get_histogram_n(), 0, 1000, 1, 10, 0 );
    pdata->histo_n_spin = gtk_spin_button_new( GTK_ADJUSTMENT(histo_n_adj), 1, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), pdata->histo_n_spin, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // Histogram size
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Bin size m" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    GtkObject *histo_m_adj = gtk_adjustment_new( _plot.get_histogram_m(), 0, 1000, 1, 10, 0 );
    pdata->histo_m_spin = gtk_spin_button_new( GTK_ADJUSTMENT(histo_m_adj), 1, 0 );
    if( type == PARTICLE_DIAG_PLOT_HISTO1D )
	gtk_widget_set_sensitive( pdata->histo_m_spin, false );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), pdata->histo_m_spin, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // Colormap interpolation style
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Bin size m" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    vbox2 = gtk_vbox_new( FALSE, 0 );
    pdata->radio_int_closest = gtk_radio_button_new_with_label_from_widget( NULL,
									     "Closest" );
    if( type == PARTICLE_DIAG_PLOT_HISTO1D )
	gtk_widget_set_sensitive( pdata->radio_int_closest, false );
    if( interpolation == INTERPOLATION_CLOSEST )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(pdata->radio_int_closest), true );
    gtk_box_pack_start( GTK_BOX(vbox2), pdata->radio_int_closest, FALSE, TRUE, 0 );
    pdata->radio_int_bilinear = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(pdata->radio_int_closest),
									      "Bilinear" );
    if( type == PARTICLE_DIAG_PLOT_HISTO1D )
	gtk_widget_set_sensitive( pdata->radio_int_bilinear, false );
    if( interpolation == INTERPOLATION_BILINEAR )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(pdata->radio_int_bilinear), true );
    gtk_box_pack_start( GTK_BOX(vbox2), pdata->radio_int_bilinear, FALSE, TRUE, 0 );
    pdata->radio_int_bicubic = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(pdata->radio_int_closest),
									     "Bicubic" );
    if( type == PARTICLE_DIAG_PLOT_HISTO1D )
	gtk_widget_set_sensitive( pdata->radio_int_bicubic, false );
    if( interpolation == INTERPOLATION_BICUBIC )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(pdata->radio_int_bicubic), true );
    gtk_box_pack_start( GTK_BOX(vbox2), pdata->radio_int_bicubic, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), vbox2, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // Dot size
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Dot size" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    GtkObject *dot_size_adj = gtk_adjustment_new( _plot.get_dot_size(), 0.1, 10.0, 0.1, 1, 0 );
    pdata->dot_size_spin = gtk_spin_button_new( GTK_ADJUSTMENT(dot_size_adj), 0.1, 1 );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), pdata->dot_size_spin, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // Ellipse fit
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Emittance ellipse fit" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    pdata->ellipse_check = gtk_check_button_new_with_label( "on/off" );
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(pdata->ellipse_check), _plot.get_emittance_ellipse() );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), pdata->ellipse_check, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // Add notebook page
    label = gtk_label_new( "Particles" );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox, label );

    return( (void *)pdata );
}


void GTKParticleDiagWindow::read_preferences( GtkWidget *notebook, void *_pdata )
{
    PreferencesData *pdata = (PreferencesData *)_pdata;

    // Plot type
    particle_diag_plot_type_e type = _plot.get_type();
    if( type != PARTICLE_DIAG_PLOT_HISTO1D ) {
	particle_diag_plot_type_e otype = type;
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(pdata->radio_plot_scatter) ) )
	    type = PARTICLE_DIAG_PLOT_SCATTER;
	else
	    type = PARTICLE_DIAG_PLOT_HISTO2D;
	if( type != otype )
	    _plot.set_type( type );
    }
    

    // Histogram size
    _plot.set_histogram_n( gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(pdata->histo_n_spin) ) );
    _plot.set_histogram_m( gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(pdata->histo_m_spin) ) );

    // Interpolation style
    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(pdata->radio_int_closest) ) )
	_plot.set_colormap_interpolation( INTERPOLATION_CLOSEST );
    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(pdata->radio_int_bilinear) ) )
	_plot.set_colormap_interpolation( INTERPOLATION_BILINEAR );
    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(pdata->radio_int_bicubic) ) )
	_plot.set_colormap_interpolation( INTERPOLATION_BICUBIC );

    // Dot size
    _plot.set_dot_size( gtk_spin_button_get_value_as_float( GTK_SPIN_BUTTON(pdata->dot_size_spin) ) );

    // Ellipse fit
    _plot.set_emittance_ellipse( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( pdata->ellipse_check) ) );

    _plot.build_plot();
    draw_and_expose();
}


void GTKParticleDiagWindow::export_data( void )
{
    GTKParticleDiagExportDialog dialog( _window, &_plot );
    dialog.run();
}


void GTKParticleDiagWindow::menuitem_export_signal( GtkToolButton *button,
						    gpointer object )
{
    GTKParticleDiagWindow *window = (GTKParticleDiagWindow *)object;
    window->export_data();
}



std::string GTKParticleDiagWindow::track_text( double x, double y )
{
    std::stringstream ss;

    particle_diag_plot_type_e type;
    trajectory_diagnostic_e diagx, diagy;
    _plot.get_plot( type, diagx, diagy );

    for( int i = 0; i < 2; i++ ) {

	double val;
	trajectory_diagnostic_e diag;
	if( i == 0 ) {
	    val = x;
	    diag = diagx;
	} else {
	    val = y;
	    diag = diagy;
	}

	switch( diag ) {
	case DIAG_X:
	    ss << "x = " << val << " m\n";
	    break;
	case DIAG_Y:
	    ss << "y = " << val << " m\n";
	    break;
	case DIAG_R:
	    ss << "r = " << val << " m\n";
	    break;
	case DIAG_Z:
	    ss << "z = " << val << " m\n";
	    break;
	case DIAG_XP:
	    ss << "x\' = " << val << " rad\n";
	    break;
	case DIAG_YP:
	    ss << "y\' = " << val << " rad\n";
	    break;
	case DIAG_RP:
	    ss << "r\' = " << val << " rad\n";
	    break;
	case DIAG_ZP:
	    ss << "z\' = " << val << " rad\n";
	    break;
	case DIAG_CURR:
	    ss << "I = " << val << " A\n";
	    break;
	default:
	    std::cout << "unknown diagnostic " << diag << "\n";
	    break;
	};
    }

    if( type == PARTICLE_DIAG_PLOT_HISTO2D ) {
	if( (diagx == DIAG_X || diagx == DIAG_Y || diagx == DIAG_R || diagx == DIAG_Z) && 
	    (diagy == DIAG_X || diagy == DIAG_Y || diagy == DIAG_R || diagy == DIAG_Z) ) {
	    // Profile plot
	    const Colormap *cmap = _plot.get_colormap();
	    double val = cmap->get_value( x, y );
	    ss << "J = " << val << " A/m2\n";
	} else {
	    // Emittance plot
	    const Colormap *cmap = _plot.get_colormap();
	    double val = cmap->get_value( x, y );
	    ss << "J = " << val << " A/(m rad)\n";
	}
    }

    
    return( ss.str() );
}





