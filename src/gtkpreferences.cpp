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













