/*! \file gtkparticlediagwindow.hpp
 *  \brief Header file for gtkparticlediagwindow.hpp
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

#ifndef GTKPARTICLEDIAGWINDOW_HPP
#define GTKPARTICLEDIAGWINDOW_HPP 1


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















