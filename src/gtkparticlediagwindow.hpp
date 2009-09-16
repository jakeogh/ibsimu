/*! \file gtkparticlediagwindow.hpp
 *  \brief Header file for interactive particle diagnostics plotter.
 */


#ifndef GTK_PARTICLE_DIAG_WINDOW_HPP
#define GTK_PARTICLE_DIAG_WINDOW_HPP 1


#include "gtkwindow.hpp"
#include "particledatabase.hpp"
#include "particlediagplot.hpp"
#include "types.hpp"


/*! \brief Interactive particle diagnostics plotter.
 */
class GTKParticleDiagWindow : public GTKWindow {

    ParticleDiagPlot _plot;

    virtual void *build_preferences( GtkWidget *notebook );
    virtual void read_preferences( GtkWidget *notebook, void *pdata );

public:

    /*! \brief Constructor for diagnostics window.
     *
     * \a style is the style of plot with 0 being scatter plot and 1
     * being colormap (histogram) plot.
     */
    GTKParticleDiagWindow( GTKPlotter *plotter, const ParticleDataBase *pdb, 
			   const Geometry *geom,
			   coordinate_axis_e axis, double level, 
			   particle_diag_plot_type_e type,
			   trajectory_diagnostic_e diagx, 
			   trajectory_diagnostic_e diagy );

    virtual ~GTKParticleDiagWindow();

};


#endif















