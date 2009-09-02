/*! \file coordmapper.hpp
 *  \brief Header file for plot coordinate mapper.
 */


#ifndef COORDMAPPER_HPP
#define COORDMAPPER_HPP 1


#include <iostream>


/*! \brief Linear 1D coordinate mapper.
 *
 *  Coordinate transformation is done like
 *
 *  x_new = xx * x + x0;
 *
 */
class Coordmapper1D {
    double _xx, _x0;   /*!< \brief Transformation coeffiecients for x-axis */

public:

    /*! \brief Default constructor for unitary transformation.
     */
    Coordmapper1D()
	: _xx(1.0), _x0(0.0) {}

    /*! \brief Constructor for fully defined transformation.
     */
    Coordmapper1D( double xx, double x0 )
	: _xx(xx), _x0(x0) {}

    /*! \brief Set transformation coefficients.
     */
    void set_transformation( double xx, double x0 ) {
	_xx = xx;
	_x0 = x0;
    }

    /*! \brief Make transformation from coordinates \a xin to
     *  coordinates \a xout.
     */
    void transform( double &xout, const double &xin ) const {
	xout = _xx * xin + _x0;
    }

    /*! \brief Make transformation from coordinates \a xin to
     *  coordinates \a xout.
     */
    void transform( double &x ) const {
	x = _xx * x + _x0;
    }

    /*! \brief Make inverse transformation for coordinate \a x.
     */
    void inv_transform( double &xout, const double &xin ) const {
	xout = (xin-_x0) / _xx;
    }

    /*! \brief Make inverse transformation for coordinate \a x.
     */
    void inv_transform( double &x ) const {
	x = (x-_x0) / _xx;
    }

    /*! \brief Debug print to stream.
     */ 
    void debug_print( std::ostream &os ) const {
	os << "**Coordmapper1D\n";
	os << "_xx = " << _xx << "\n";
	os << "_x0 = " << _x0 << "\n";
    }
};



/*! \brief Linear-linear 2D coordinate mapper.
 *
 */
class Coordmapper {
    Coordmapper1D _cmx;  /*!< \brief Transformation for x-axis */
    Coordmapper1D _cmy;  /*!< \brief Transformation for y-axis */

public:

    /*! \brief Default constructor for unitary transformation.
     */
    Coordmapper() {}

    /*! \brief Constructor for fully defined transformation using 1D
     *  coordinate mappers.
     */
    Coordmapper( Coordmapper1D cmx, Coordmapper1D cmy )
	: _cmx(cmx), _cmy(cmy) {}

    /*! \brief Constructor for fully defined transformation.
     */
    Coordmapper( double xx, double x0, double yy, double y0 )
	: _cmx(xx,x0), _cmy(yy,y0) {}

    /*! \brief Set transformation matrix coefficients.
     */
    void set_transformation( double xx, double x0, double yy, double y0 ) {
	_cmx.set_transformation( xx, x0 );
	_cmy.set_transformation( yy, y0 );
    }

    /*! \brief Make transformation for coordinates \a x, \a y.
     */
    void transform( double &x, double &y ) const {
	_cmx.transform( x );
	_cmy.transform( y );
    }

    /*! \brief Make transformation from coordinates \a xin to
     *  coordinates \a xout.
     */
    void transform( double *xout, const double *xin ) const {
	_cmx.transform( xout[0], xin[0] );
	_cmy.transform( xout[1], xin[1] );
    }

    /*! \brief Make inverse transformation for coordinates \a x, \a y.
     */
    void inv_transform( double &x, double &y ) const {
	_cmx.inv_transform( x );
	_cmy.inv_transform( y );
    }

    /*! \brief Make inverse transformation from coordinates \a xin to
     *  coordinates \a xout.
     */
    void inv_transform( double *xout, const double *xin ) const {
	_cmx.inv_transform( xout[0], xin[0] );
	_cmy.inv_transform( xout[1], xin[1] );
    }

};


#endif















