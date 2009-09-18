/*! \file gtkpreferences.cpp
 *  \brief Source code for gtkpreferences.cpp
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

#include "gtkpreferences.hpp"
#include "gtkwindow.hpp"


GTKPreferences::GTKPreferences( GTKWindow *gtkwindow, GtkWidget *window, Frame *frame )
    : _gtkwindow(gtkwindow), _window(window), _frame(frame)
{

}


GTKPreferences::~GTKPreferences()
{

}


void GTKPreferences::run( void )
{
    GtkWidget *dialog = gtk_dialog_new_with_buttons( "Plot preferences",
						     GTK_WINDOW(_window),
						     (GtkDialogFlags)(GTK_DIALOG_MODAL | 
								      GTK_DIALOG_DESTROY_WITH_PARENT), 
						     GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
						     GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
						     NULL );
    GtkWidget *vbox = GTK_DIALOG(dialog)->vbox;

    // ****************************************************************************

    GtkWidget *vbox2 = gtk_vbox_new( FALSE, 0 );

    // Fontsize
    GtkWidget *hbox = gtk_hbox_new( TRUE, 30 );
    GtkWidget *label = gtk_label_new( "Fontsize" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    GtkWidget *fontsize_entry = gtk_entry_new();
    std::string s = to_string( _frame->get_font_size() );
    gtk_entry_set_text( GTK_ENTRY(fontsize_entry), s.c_str() );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), fontsize_entry, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox2), hbox, FALSE, TRUE, 0 );

    // Range xmin
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Range xmin" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    GtkWidget *rxmin_entry = gtk_entry_new();
    double min, max;
    _frame->get_ranges( PLOT_AXIS_X1, min, max );
    s = to_string( min );
    gtk_entry_set_text( GTK_ENTRY(rxmin_entry), s.c_str() );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), rxmin_entry, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox2), hbox, FALSE, TRUE, 0 );
    
    // Range xmax
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Range xmax" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    GtkWidget *rxmax_entry = gtk_entry_new();
    _frame->get_ranges( PLOT_AXIS_X1, min, max );
    s = to_string( max );
    gtk_entry_set_text( GTK_ENTRY(rxmax_entry), s.c_str() );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), rxmax_entry, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox2), hbox, FALSE, TRUE, 0 );

    // Range ymin
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Range ymin" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    GtkWidget *rymin_entry = gtk_entry_new();
    _frame->get_ranges( PLOT_AXIS_Y1, min, max );
    s = to_string( min );
    gtk_entry_set_text( GTK_ENTRY(rymin_entry), s.c_str() );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), rymin_entry, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox2), hbox, FALSE, TRUE, 0 );
    
    // Range ymax
    hbox = gtk_hbox_new( TRUE, 30 );
    label = gtk_label_new( "Range ymax" );
    gtk_misc_set_alignment( GTK_MISC(label), 0, 0.5 );
    GtkWidget *rymax_entry = gtk_entry_new();
    _frame->get_ranges( PLOT_AXIS_Y1, min, max );
    s = to_string( max );
    gtk_entry_set_text( GTK_ENTRY(rymax_entry), s.c_str() );
    gtk_box_pack_start( GTK_BOX(hbox), label, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), rymax_entry, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox2), hbox, FALSE, TRUE, 0 );

    // Notebook, page 1
    label = gtk_label_new( "Frame" );
    GtkWidget *notebook = gtk_notebook_new();
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );

    // Notebook, additional pages
    void *pdata = _gtkwindow->build_preferences( notebook );

    // Pack notebook
    gtk_box_pack_start( GTK_BOX(vbox), notebook, FALSE, TRUE, 0 );

    gtk_widget_show_all( dialog );
    if( gtk_dialog_run( GTK_DIALOG(dialog) ) == GTK_RESPONSE_ACCEPT ) {

	// Fontsize
	double fontsize = atof( gtk_entry_get_text( GTK_ENTRY(fontsize_entry) ) );
	_frame->set_font_size( fontsize );

	// Ranges
	double rxmin = atof( gtk_entry_get_text( GTK_ENTRY(rxmin_entry) ) );
	double rxmax = atof( gtk_entry_get_text( GTK_ENTRY(rxmax_entry) ) );
	double rymin = atof( gtk_entry_get_text( GTK_ENTRY(rymin_entry) ) );
	double rymax = atof( gtk_entry_get_text( GTK_ENTRY(rymax_entry) ) );
	_frame->set_ranges( PLOT_AXIS_X1, rxmin, rxmax );
	_frame->set_ranges( PLOT_AXIS_Y1, rymin, rymax );

	// Read additional pages of notebook
	_gtkwindow->read_preferences( notebook, pdata );

	// Refresh frame
	_gtkwindow->draw_and_expose();
    }

    gtk_widget_destroy( dialog );
}













