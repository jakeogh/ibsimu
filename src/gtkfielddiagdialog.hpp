#include <gtk/gtk.h>

#include "gtkplotter.hpp"


#ifndef GTK_FIELD_DIAG_DIALOG_HPP
#define GTK_FIELD_DIAG_DIALOG_HPP 1


class GTKFieldDiagDialog
{

    GtkWidget      *_window;
    GTKPlotter     *_plotter;
    const Geometry *_geom;

    double          _x1[3];
    double          _x2[3];

public:

    GTKFieldDiagDialog( GtkWidget *window, GTKPlotter *plotter, double x1[3], double x2[3] );

    ~GTKFieldDiagDialog();

    void run( void );
};


#endif













