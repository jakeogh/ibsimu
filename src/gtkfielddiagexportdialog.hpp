#include <gtk/gtk.h>

#include "gtkplotter.hpp"


#ifndef GTK_FIELD_DIAG_EXPORT_DIALOG_HPP
#define GTK_FIELD_DIAG_EXPORT_DIALOG_HPP 1


/*! \brief Interactive dialog for exporting field diagnostic data.
 */
class GTKFieldDiagExportDialog
{
    GtkWidget           *_window;
    const FieldDiagPlot *_plot;

public:

    GTKFieldDiagExportDialog( GtkWidget *window, const FieldDiagPlot *plot );

    ~GTKFieldDiagExportDialog();

    void run( void );
};


#endif













