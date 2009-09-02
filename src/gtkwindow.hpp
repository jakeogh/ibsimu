#ifndef GTK_WINDOW_HPP
#define GTK_WINDOW_HPP 1


#include <cairo.h>
#include <gtk/gtk.h>

#include "gtkplotter.hpp"
#include "frame.hpp"
#include "graph3d.hpp"



class GTKWindow {

    size_t            _width;
    size_t            _height;

    GtkWidget        *_trackwindow;
    GtkWidget        *_tracklabel;

    cairo_t          *_cairo;
    cairo_surface_t  *_surface;

    int               _start[2];
    int               _end[2];
    int               _tool;

    void hardcopy( void );
    void move( int action, double x, double y );
    void track( int action, double x, double y );
    void zoom_out( double x, double y );
    void zoom_in( double x, double y );
    void zoom_window( int action, double x, double y );

    void darea_motion( GdkEventMotion *event );
    void darea_enter( GdkEventCrossing *event );
    void darea_leave( GdkEventCrossing *event );
    void darea_button( GdkEventButton *event );
    void configure( void );
    void derea_button( GdkEventButton *event );
    void delete_window( void );
    void menuitem_tool_change( GtkToolButton *button );
    void menuitem_preferences( GtkMenuItem *menuitem );

    static gboolean window_delete_signal( GtkWidget *widget, 
					  GdkEventExpose *event, 
					  gpointer object );

    static void menuitem_quit_signal( GtkMenuItem *menuitem,
				      gpointer object );
    static void menuitem_tool_change_signal( GtkToolButton *button,
					     gpointer object );
    static void menuitem_hardcopy_signal( GtkToolButton *button,
					  gpointer object );
    static void menuitem_zoom_fit_signal( GtkToolButton *button,
					  gpointer object );
    static void menuitem_preferences_signal( GtkMenuItem *menuitem,
					     gpointer object );

    static gboolean darea_configure_signal( GtkWidget *widget, 
					    GdkEventConfigure *event, 
					    gpointer object );
    static gboolean darea_expose_signal( GtkWidget *widget, 
					 GdkEventExpose *event, 
					 gpointer object );
    static gboolean darea_button_signal( GtkWidget *widget, 
					 GdkEventButton *event,
					 gpointer object );
    static gboolean darea_motion_signal( GtkWidget *widget, 
					 GdkEventMotion *event,
					 gpointer object );
    static gboolean darea_enter_signal( GtkWidget *widget, 
					GdkEventCrossing *event,
					gpointer object );
    static gboolean darea_leave_signal( GtkWidget *widget, 
					GdkEventCrossing *event,
					gpointer object );

protected:

    GTKPlotter       *_plotter;

    Frame             _frame;
    GtkWidget        *_window;
    GtkWidget        *_darea;
    GtkWidget        *_menubar;
    GtkWidget        *_toolbar;
    GtkWidget        *_statusbar;

    GtkToolItem      *_radioitem;

    GTKWindow( GTKPlotter *plotter );

    virtual void zoom_fit( void );
    virtual std::string track_text( double x, double y );
    void expose( int x, int y, int width, int height );

public:

    virtual ~GTKWindow();

    virtual void *build_preferences( GtkWidget *notebook );
    virtual void read_preferences( GtkWidget *notebook, void *pdata );

    void draw_and_expose( void );
    void show( void );
};





#endif













