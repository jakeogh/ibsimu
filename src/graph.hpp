/*! \file graph.hpp
 *  \brief Header file for plot drawable abstract base class
 */


#ifndef DRAWABLE_HPP
#define DRAWABLE_HPP 1


#include <cairo.h>
#include "coordmapper.hpp"


/*! \brief Abstract base class for drawable plots.
 */
class Graph {

public:

    /*! \brief Virtual destructor.
     */
    virtual ~Graph() {}

    /*! \brief Plot graph with cairo.
     *
     *  Plot the graph using \a cairo and coordinate mapper \a
     *  cm. The visible range of plot is given in array \a range in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] ) = 0;

    /*! \brief Get bounding box of graph.
     *
     *  Returns the bounding box of the graph in array \a bbox in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void get_bbox( double bbox[4] ) = 0;
};


#endif













