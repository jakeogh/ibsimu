/*! \file graph3d.hpp
 *  \brief Header file for plot drawable3d abstract base class
 */


#ifndef DRAWABLE3D_HPP
#define DRAWABLE3D_HPP 1


#include "graph.hpp"


enum view_e {
    VIEW_XY = 0,
    VIEW_XZ,
    VIEW_YZ
};


/*! \brief Abstract base class for drawable3d plots.
 */
class Graph3D : public Graph {

protected:

    view_e     _view;
    int        _vb[3];
    int        _level;

public:

    /*! \brief Constructor.
     */
    Graph3D() {
	_view  = VIEW_XY;
	_vb[0] = 0;
	_vb[1] = 1;
	_vb[2] = 2;
	_level = 0;
    }

    /*! \brief Virtual destructor.
     */
    virtual ~Graph3D() {}

    /*! \brief Plot drawable with cairo.
     *
     *  Plot the drawable using \a cairo and coordinate mapper \a
     *  cm. The visible range of plot is given in array \a range in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] ) = 0;

    /*! \brief Get bounding box of drawable.
     *
     *  Returns the bounding box of the drawable in array \a bbox in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void get_bbox( double bbox[4] ) = 0;

    /*! \brief Set the view of 3D drawable.
     *
     *  Sets view direction to \a view and the view level to \a level.
     */
    void set_view( view_e view, int level ) {
        if( view == VIEW_XY ) {
	    _vb[0] = 0;
	    _vb[1] = 1;
	    _vb[2] = 2;
	} else if( view == VIEW_XZ ) {
	    _vb[0] = 0;
	    _vb[1] = 2;
	    _vb[2] = 1;
	} else if( view == VIEW_YZ ) {
	    _vb[0] = 1;
	    _vb[1] = 2;
	    _vb[2] = 0;
	}
	_view = view;
	_level = level;
    }
};


#endif













