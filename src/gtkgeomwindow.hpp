#ifndef GTKGEOMWINDOW_HPP
#define GTKGEOMWINDOW_HPP 1


#include <vector>

#include "gtkwindow.hpp"
#include "geomplot.hpp"
#include "vectorfield.hpp"



/*! \brief Interactive geometry plotter.
 */
class GTKGeomWindow : public GTKWindow {

    GeomPlot                _geomplot;

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













