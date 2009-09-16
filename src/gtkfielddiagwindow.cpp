#include <limits>
#include "gtkfielddiagwindow.hpp"
#include "gtkfielddiagexportdialog.hpp"


GTKFieldDiagWindow::GTKFieldDiagWindow( GTKPlotter *plotter, const Geometry *geom, size_t N, 
					const Vec3D &x1, const Vec3D &x2,
					const field_diag_type_e diag[2], 
					const field_loc_type_e loc[2] )
    : GTKWindow(plotter), _plot(&_frame,geom)
{
    _plot.set_epot( plotter->get_epot() );
    _plot.set_efield( plotter->get_efield() );
    _plot.set_scharge( plotter->get_scharge() );
    _plot.set_bfield( plotter->get_bfield() );

    _plot.set_coordinates( N, x1, x2 );
    _plot.set_diagnostic( diag, loc );

    // Set window title
    gtk_window_set_title( GTK_WINDOW(_window), "Field diagnostics" );

    // Add export menu item
    GtkWidget *item_export = gtk_menu_item_new_with_mnemonic( "_Export" );
    gtk_menu_shell_prepend( GTK_MENU_SHELL(_menu_file), item_export );
    g_signal_connect( G_OBJECT(item_export), "activate",
		      G_CALLBACK(menuitem_export_signal),
		      (gpointer)this );

    _plot.build_plot();
    _frame.get_ranges( PLOT_AXIS_X1, _x1min, _x1max );
    _frame.get_ranges( PLOT_AXIS_X1, _x2min, _x2max );
    show();
}


GTKFieldDiagWindow::~GTKFieldDiagWindow()
{
}


void GTKFieldDiagWindow::export_data( void )
{
    GTKFieldDiagExportDialog dialog( _window, &_plot );
    dialog.run();
}


void GTKFieldDiagWindow::menuitem_export_signal( GtkToolButton *button,
						 gpointer object )
{
    GTKFieldDiagWindow *window = (GTKFieldDiagWindow *)object;
    window->export_data();
}


void GTKFieldDiagWindow::zoom_fit( void )
{
    double min = -std::numeric_limits<double>::infinity();
    double max = std::numeric_limits<double>::infinity();
    _frame.set_ranges( PLOT_AXIS_X1, _x1min, _x1max );
    _frame.set_ranges( PLOT_AXIS_Y1, min, max );
    _frame.set_ranges( PLOT_AXIS_X2, _x2min, _x2max );
    _frame.set_ranges( PLOT_AXIS_Y2, min, max );
    _frame.ruler_autorange_enable( PLOT_AXIS_X1, false, false );
    _frame.ruler_autorange_enable( PLOT_AXIS_Y1, false, false );
    _frame.ruler_autorange_enable( PLOT_AXIS_X2, false, false );
    _frame.ruler_autorange_enable( PLOT_AXIS_Y2, false, false );

    draw_and_expose();
}


void *GTKFieldDiagWindow::build_preferences( GtkWidget *notebook )
{
    return( NULL );
}


void GTKFieldDiagWindow::read_preferences( GtkWidget *notebook, void *pdata )
{

}













