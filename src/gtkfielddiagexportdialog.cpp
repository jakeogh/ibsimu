#include "gtkfielddiagexportdialog.hpp"


GTKFieldDiagExportDialog::GTKFieldDiagExportDialog( GtkWidget *window, const FieldDiagPlot *plot )
    : _window(window), _plot(plot)
{

}


GTKFieldDiagExportDialog::~GTKFieldDiagExportDialog()
{

}


void GTKFieldDiagExportDialog::run( void )
{
   GtkWidget *dialog = gtk_file_chooser_dialog_new( "Export field data",
						    GTK_WINDOW(_window),
						    GTK_FILE_CHOOSER_ACTION_SAVE,
						    GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
						    GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
						    NULL );

   gtk_file_chooser_set_current_name( GTK_FILE_CHOOSER(dialog), "field.txt" );
   gtk_file_chooser_set_show_hidden( GTK_FILE_CHOOSER(dialog), TRUE );
   gtk_file_chooser_set_local_only(GTK_FILE_CHOOSER(dialog), TRUE );
   
   gtk_widget_show_all( dialog );
   if( gtk_dialog_run( GTK_DIALOG(dialog) ) == GTK_RESPONSE_ACCEPT ) {

       char *filename = gtk_file_chooser_get_filename( GTK_FILE_CHOOSER(dialog) );
       std::cout << filename << "\n";

       // Write output to filename
       _plot->export_data( filename );

       g_free( filename );
   }

   gtk_widget_destroy( dialog );
}
