#include <gtk/gtk.h>

#include "gtkplotter.hpp"


#ifndef GTK_PARTICLE_DIAG_DIALOG_HPP
#define GTK_PARTICLE_DIAG_DIALOG_HPP 1


class GTKParticleDiagDialog
{

    GtkWidget      *_window;
    GTKPlotter     *_plotter;
    const Geometry *_geom;

    int             _plane;
    double          _val;


    GtkWidget      *_radio_plane_x;
    GtkWidget      *_radio_plane_y;
    GtkWidget      *_radio_plane_z;

    GtkWidget      *_radio_emit_xx;
    GtkWidget      *_radio_emit_yy;
    GtkWidget      *_radio_emit_zz;
    GtkWidget      *_radio_prof_yz; // Or prof_x in case of 2d
    GtkWidget      *_radio_prof_xz; // Or prof_y in case of 2d
    GtkWidget      *_radio_prof_xy;

    GtkWidget      *_radio_plot_scatter;
    GtkWidget      *_radio_plot_colormap;

    void profile_toggled2( GtkToggleButton *togglebutton );
    static void profile_toggled( GtkToggleButton *togglebutton,
				 gpointer         user_data );

    void plane_activated( void );
    static void plane_toggled( GtkToggleButton *togglebutton,
			       gpointer         user_data );


public:

    GTKParticleDiagDialog( GtkWidget *window, GTKPlotter *plotter, int plane, double val );

    ~GTKParticleDiagDialog();

    void run( void );
};


#endif













