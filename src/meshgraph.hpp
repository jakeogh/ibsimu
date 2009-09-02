/*! \file meshgraph.hpp
 *  \brief Header file for mesh plotter.
 */


#ifndef MESHPLOT_HPP
#define MESHPLOT_HPP 1


#include "graph3d.hpp"
#include "geometry.hpp"


/*! \brief Class for drawing mesh
 *
 *  Implementation of %Graph3D.
 */
class MeshGraph : public Graph3D {

    const Geometry         &_g;               /*!< \brief Reference to simulation geometry. */

public:

    /*! \brief Constructor for mesh plotter.
     */
    MeshGraph( const Geometry &g );

    /*! \brief Destructor.
     */
    virtual ~MeshGraph();

    /*! \brief Plot graph with cairo.
     *
     *  Plot the graph using \a cairo and coordinate mapper \a
     *  cm. The visible range of plot is given in array \a range in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] );

    /*! \brief Get bounding box of graph.
     *
     *  Returns the bounding box of the graph in array \a bbox in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void get_bbox( double bbox[4] );
};


#endif













