/*! \file lineclip.hpp
 *  \brief Header file for line clipper.
 */


#ifndef LINECLIP_HPP
#define LINECLIP_HPP 1

#include <cairo.h>


/*! \brief Line clipper.
 *
 *  Cairo graphics coordinates are internally handled using fixed
 *  point algebra for speed. This causes problems in clipping
 *  algorithm when there are large scale differences in
 *  coordinates. This class is provided for the user to overcomme this
 *  problem by using floating point algebra for line clipping.
 */
class LineClip {
    cairo_t  *p_dc;          /* Cairo context */
    double    clip[4];       /* xmin, ymin, xmax, ymax  */
    
    double    first[2];      /* Last user given moveto point (first point of path) */
    
    double    last[2];       /* last user given point, nan if not available */
    int       last_outcode;  /* outcode of last user given point */
    int       last_op;       /* last operation 0=lineto, 1=moveto, 2=no op */
    
    double    drawn[2];      /* last drawn point */
    int       drawn_outcode; /* outcode of last drawn point */
    
    int       coord_alloc;   /* Allocated size of coordinate database */
    double   *coord;         /* Buffer for storing line coordinates in curve_to */

    int outcode( double x, double y );
    int exit_outcode( double x, double y );
    void get_point( double *coords, double t,
		    double x0, double y0,
		    double x1, double y1,
		    double x2, double y2,
		    double x3, double y3 );

public:

    LineClip( cairo_t *cairo );
    ~LineClip();

    void set( double xmin, double ymin, double xmax, double ymax );
    void reset();

    void move_to( double x, double y );
    void line_to( double x, double y );
    void curve_to( double x1, double y1,
		   double x2, double y2,
		   double x3, double y3 );
    void close_path();
    void fill();
};


#endif













