/*! \file histogram.hpp
 *  \brief %Histogram data handling for 1D and 2D
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * taneli.kalvas@jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */

#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP 1


#include <vector>


/*! \brief Base histogram class.
 */
class Histogram
{

public:

    /*! \brief Destructor.
     */
    virtual ~Histogram() {}

};


/*! \brief Class for 1D histogram type representation of data.
 */
class Histogram1D : public Histogram
{
    int                 _n;         /*!< \brief Number of bins. */
    double              _range[2];  /*!< \brief Ranges: min, max. */
    double              _step;      /*!< \brief Step size. */
    std::vector<double> _data;      /*!< \brief Data of histogram. */

public:

    /*! \brief Constructor for \a n bin histogram with \a ranges.
     */
    Histogram1D( size_t n, const double range[2] );

    /*! \brief Constructor for \a n bin histogram from scatter data with even weights.
     */
    Histogram1D( size_t n, const std::vector<double> &xdata );

    /*! \brief Constructor for \a n bin histogram from scatter data with weights wrom \a wdata.
     */
    Histogram1D( size_t n, const std::vector<double> &xdata, const std::vector<double> &wdata );

    /*! \brief Destructor.
     */
    virtual ~Histogram1D();

    /*! \brief Return the number of bins.
     */
    size_t n( void ) const { return( _n ); }

    /*! \brief Return the step size.
     */
    double step( void ) const { return( _step ); }

    /*! \brief Return the coordinate on bin \a i.
     */
    double coord( size_t i ) const { return( _range[0] + i*(_range[1]-_range[0]) / (_n-1.0) ); }

    /*! \brief Accumulate \a weight on bin \a i.
     *
     *  Not a safe function. Input not checked.
     */
    void accumulate( size_t i, double weight ) {
	_data[i] += weight;
    }

    /*! \brief Accumulate \a weight on bins around \a x linearly.
     *
     *  Accumulation is done on two neighbouring bins around point \a
     *  x. The distribution of weight is done using inverse linear
     *  interpolation.
     *
     *  This is a safe function. Accumulation outside histogram range
     *  is discarded.
     */
    void accumulate_linear( double x, double weight );

    /*! \brief Return data range.
     */
    void get_range( double range[2] ) const { 
	range[0] =  _range[0];
	range[1] =  _range[1];
    }
    
    /*! \brief Return bin range.
     *
     *  Returns minimum and maximum values on any bin in histogram.
     */
    void get_bin_range( double &min, double &max ) const;
    
    /*! \brief Return a reference to the histogram data.
     */
    std::vector<double> &get_data( void ) { return( _data ); }

    /*! \brief Return a reference to the histogram data.
     */
    const std::vector<double> &get_data( void ) const { return( _data ); }

    /*! \brief Return a const reference to the weight on bin \a i.
     */
    const double &operator()( size_t i ) const {
	return( _data[i] );
    }

    /*! \brief Return a reference to the weight on bin \a i.
     */
    double &operator()( size_t i ) {
	return( _data[i] );
    }
};


/*! \brief Class for 2d histogram type representation of data.
 */
class Histogram2D : public Histogram
{
    int                 _n;         /*!< \brief Number of bins along first axis. */
    int                 _m;         /*!< \brief Number of bins along second axis. */
    double              _range[4];  /*!< \brief Ranges: Amin, Bmin, Amax, Bmax. */
    double              _nstep;     /*!< \brief Step size along first axis. */
    double              _mstep;     /*!< \brief Step size along first axis. */
    std::vector<double> _data;      /*!< \brief Data of histogram. */

public:

    /*! \brief Constructor for \a n x \a m histogram with \a ranges.
     */
    Histogram2D( size_t n, size_t m, const double range[4] );

    /*! \brief Constructor for \a n x \a m histogram from scatter xy-data with even weights.
     */
    Histogram2D( size_t n, size_t m, 
		 const std::vector<double> &xdata,
		 const std::vector<double> &ydata );

    /*! \brief Constructor for \a n x \a m histogram from scatter xy-data with weights from \a wdata.
     */
    Histogram2D( size_t n, size_t m, 
		 const std::vector<double> &xdata,
		 const std::vector<double> &ydata,
		 const std::vector<double> &wdata );

    /*! \brief Destructor.
     */
    virtual ~Histogram2D();

    /*! \brief Return the number of bins along the first axis.
     */
    size_t n( void ) const { return( _n ); }

    /*! \brief Return the number of bins along the second axis.
     */
    size_t m( void ) const { return( _m ); }

    /*! \brief Return the step size along along the first axis.
     */
    double nstep( void ) const { return( _nstep ); }

    /*! \brief Return the step size along along the second axis.
     */
    double mstep( void ) const { return( _mstep ); }

    /*! \brief Return the coordinate along the first axis on bin \a i.
     */
    double icoord( size_t i ) const { return( _range[0] + i*(_range[2]-_range[0]) / (_n-1.0) ); }

    /*! \brief Return the coordinate along the second axis on bin \a j.
     */
    double jcoord( size_t j ) const { return( _range[1] + j*(_range[3]-_range[1]) / (_m-1.0) ); }

    /*! \brief Accumulate \a weight on bin \a (i,j).
     *
     *  Not a safe function. Input not checked.
     */
    void accumulate( size_t i, size_t j, double weight ) {
	_data[i+j*_n] += weight;
    }

    /*! \brief Accumulate \a weight on bins around \a (x,y) linearly.
     *
     *  Accumulation is done on four neighbouring bins around point \a
     *  (x,y). The distribution of weight is done using inverse bilinear
     *  interpolation.
     *
     *  This is a safe function. Accumulation outside histogram range
     *  is discarded.
     */
    void accumulate_linear( double x, double y, double weight );

    /*! \brief Return data range.
     */
    void get_range( double range[4] ) const { 
	range[0] =  _range[0];
	range[1] =  _range[1];
	range[2] =  _range[2];
	range[3] =  _range[3];
    }
    
    /*! \brief Return bin range.
     *
     *  Returns minimum and maximum values on any bin in histogram.
     */
    void get_bin_range( double &min, double &max ) const;
    
    /*! \brief Return a reference to the histogram data.
     */
    std::vector<double> &get_data( void ) { return( _data ); }

    /*! \brief Return a reference to the histogram data.
     */
    const std::vector<double> &get_data( void ) const { return( _data ); }

    /*! \brief Return a const reference to the weight on bin \a (i,j).
     */
    const double &operator()( size_t i, size_t j ) const {
	return( _data[i+j*_n] );
    }

    /*! \brief Return a reference to the weight on bin \a (i,j).
     */
    double &operator()( size_t i, size_t j ) {
	return( _data[i+j*_n] );
    }
};


#endif


















