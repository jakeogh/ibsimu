/*! \file fieldgraph.hpp
 *  \brief Header file for mesh plotter.
 */


#ifndef FIELDPLOT_HPP
#define FIELDPLOT_HPP 1


#include <vector>
#include "graph3d.hpp"
#include "scalarfield.hpp"
#include "colormap.hpp"
#include "geometry.hpp"


/*! \brief Class for drawing fields with colormap
 *
 *  Implementation of %Graph3D.
 */
class FieldGraph : public Graph3D {

    const ScalarField      *_scalarfield;     /*!< \brief Scalarfield for plotting. */
    Colormap               *_colormap;        /*!< \brief Colormap for field plot. */

    view_e                  _oview;
    double                  _olevel;

    bool                    _enabled;         /*!< \brief Is mesh plotting enabled */

    void build_scalarfield_plot( void );

public:

    /*! \brief Constructor for plotting ScalarField.
     */
    FieldGraph( const ScalarField *field );

    /*! \brief Destructor.
     */
    virtual ~FieldGraph();

    /*! \brief Enable/disable plot.
     */
    void enable( bool enable );

    /*! \brief Plot drawable with cairo.

     *  Plot drawable using \a cairo and coordinate mapper \a cm. The
     *  visible range of plot is given in array \a range in order
     *  xmin, ymin, xmax, ymax.
     */
    virtual void plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] );

    /*! \brief Get bounding box of drawable.
     *
     *  Returns the bounding box of the drawable in array \a bbox in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void get_bbox( double bbox[4] );
};


#endif













