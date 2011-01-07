/*! \file gtkgeomwindow.hpp
 *  \brief %Geometry view window
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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
 * taneli.kalvas@jyu.fi.
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

#ifndef GTKGEOMWINDOW_HPP
#define GTKGEOMWINDOW_HPP 1


#include <vector>

#include "gtkwindow.hpp"
#include "geomplot.hpp"
#include "vectorfield.hpp"



/*! \brief Interactive geometry plotter.
 */
class GTKGeomWindow : public GTKWindow {

    GeomPlot                 _geomplot;

    const Geometry          *_geom;
    const ScalarField       *_epot;
    const ScalarField       *_scharge;
    const VectorField       *_bfield;
    const ParticleDataBase  *_pdb;

    int                      _tool;
    int                      _start[2];
    int                      _end[2];

    GtkWidget               *_spinbutton;
    GtkWidget               *_combobox;


    void update_view();

    virtual void zoom_fit( void );
    virtual std::string track_text( double x, double y );
 
    virtual void *build_preferences( GtkWidget *notebook );
    virtual void read_preferences( GtkWidget *notebook, void *pdata );

    void combobox( GtkComboBox *combobox );
    void spinbutton( GtkSpinButton *spinbutton );
    void menuitem_tool_change( GtkToolButton *button );
    void field_diag( int action, double x, double y );
    void particle_diag( int action, double x, double y );
    void darea_motion2( GdkEventMotion *event );
    void darea_button2( GdkEventButton *event );

    static void combobox_signal( GtkComboBox *combobox,
				 gpointer object );
    static void spinbutton_signal( GtkSpinButton *spinbutton,
				   gpointer object );
    static void menuitem_tool_change_signal( GtkToolButton *button,
					     gpointer object );
    static gboolean darea_motion_signal2( GtkWidget *widget, 
					  GdkEventMotion *event,
					  gpointer object );
    static gboolean darea_button_signal2( GtkWidget *widget, 
					  GdkEventButton *event,
					  gpointer object );


public:
    
    GTKGeomWindow( class GTKPlotter *plotter,
		   const Geometry *geom,
		   const ScalarField *epot,
		   const ScalarField *scharge,
		   const VectorField *bfield,
		   const ParticleDataBase *pdb );
    
    ~GTKGeomWindow();
};






#endif


















