#ifndef GTK_FIELD_DIAG_WINDOW_HPP
#define GTK_FIELD_DIAG_WINDOW_HPP 1


#include "gtkwindow.hpp"
#include "fielddiagplot.hpp"


class GTKFieldDiagWindow : public GTKWindow {

    FieldDiagPlot              _plot;

    double                     _x1min;
    double                     _x1max;
    double                     _x2min;
    double                     _x2max;

    virtual void zoom_fit( void );
    virtual void *build_preferences( GtkWidget *notebook );
    virtual void read_preferences( GtkWidget *notebook, void *pdata );

public:

    /*! \brief Make new field diagnostic window.
     *
     *  Field is diagnosed on a line from \a x1 to \a x2. On the first
     *  y-axis, \a g1type field is plotted. On the second \a g2type
     *  field is plotted. The data displayed on first x-axis is
     *  specified by \a dist1type and the second x-axis on \a
     *  dist2type.
     */
    GTKFieldDiagWindow( GTKPlotter *plotter, const Geometry *geom, size_t N, 
			const Vec3D &x1, const Vec3D &x2, 
			const field_diag_type_e diag[2], const field_loc_type_e loc[2] );

    virtual ~GTKFieldDiagWindow();

};


#endif















