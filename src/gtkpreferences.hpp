/*! \file gtkpreferences.hpp
 *  \brief Header file for simulation preferences.
 */


#ifndef GTK_PREFERENCES_HPP
#define GTK_PREFERENCES_HPP 1


#include <gtk/gtk.h>
#include <png.h>
#include "frame.hpp"



/*! \brief Common features of preferences settings for all gtk user
 *  interface based interactive plotters.
 */
class GTKPreferences
{
    class GTKWindow *_gtkwindow;
    GtkWidget       *_window;
    Frame           *_frame;

public:

    GTKPreferences( class GTKWindow *gtkwindow, GtkWidget *window, Frame *frame );

    ~GTKPreferences();

    void run( void );

};


#endif













