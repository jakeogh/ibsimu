#ifndef FIELDDIAGPLOTTER_HPP
#define FIELDDIAGPLOTTER_HPP 1


#include "geometry.hpp"

#include "fielddiagplot.hpp"
#include "plotter.hpp"



/*! \brief Non-interactive fielddiagetry plotter.
 *
 *  This class ties together Plotter, which provides basic graphics
 *  features including frame control, and FielddiagPlot, which makes the
 *  fielddiagetry plot.
 */
class FieldDiagPlotter : public Plotter, public FieldDiagPlot {

    virtual void build_plot( void );

public:

    /*! \brief Constructor for field diagnostic plotter.
     */
    FieldDiagPlotter( const Geometry *geom );

    /*! \brief Destructor for field diagnostic plotter.
     */
    ~FieldDiagPlotter();


};



#endif

