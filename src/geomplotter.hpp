#ifndef GEOMPLOTTER_HPP
#define GEOMPLOTTER_HPP 1


#include "geometry.hpp"

#include "geomplot.hpp"
#include "plotter.hpp"



/*! \brief Non-interactive geometry plotter.
 *
 *  This class ties together Plotter, which provides basic graphics
 *  features including frame control, and GeomPlot, which makes the
 *  geometry plot.
 */
class GeomPlotter : public Plotter, public GeomPlot {

    virtual void build_plot( void );

public:

    /*! \brief Constructor for geometry plotter.
     */
    GeomPlotter( const Geometry *geom );

    /*! \brief Destructor for geometry plotter.
     */
    ~GeomPlotter();


};



#endif
