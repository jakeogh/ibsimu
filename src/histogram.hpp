/*! \file histogram.hpp
 *  \brief Header file for histograms.
 */


#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP 1


#include <vector>


/*! \brief Class for 1D histogram type representation of data.
 */
class Histogram
{
    int                 _n;         /*!< \brief Number of bins. */
    double              _range[2];  /*!< \brief Ranges: min, max. */
    double              _step;      /*!< \brief Step size. */
    std::vector<double> _data;      /*!< \brief Data of histogram. */

public:

    /*! \brief Constructor for \a n bin histogram with \a ranges.
     */
    Histogram( size_t n, const double range[2] );

    /*! \brief Constructor for \a n bin histogram from scatter data with even weights.
     */
    Histogram( size_t n, const std::vector<double> &xdata );

    /*! \brief Constructor for \a n bin histogram from scatter data with weights wrom \a wdata.
     */
    Histogram( size_t n, const std::vector<double> &xdata, const std::vector<double> &wdata );

    /*! \brief Destructor.
     */
    ~Histogram();

    /*! \brief Return the number of bins.
     */
    size_t n( void ) { return( _n ); }

    /*! \brief Return the step size.
     */
    double step( void ) { return( _step ); }

    /*! \brief Return the coordinate on bin \a i.
     */
    double coord( size_t i ) { return( _range[0] + i*(_range[1]-_range[0]) / (_n-1.0) ); }

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
    void get_range( double range[2] ) { 
	range[0] =  _range[0];
	range[1] =  _range[1];
    }
    
    /*! \brief Return bin range.
     *
     *  Returns minimum and maximum values on any bin in histogram.
     */
    void get_bin_range( double &min, double &max );
    
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
class Histogram2D
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
    ~Histogram2D();

    /*! \brief Return the number of bins along the first axis.
     */
    size_t n( void ) { return( _n ); }

    /*! \brief Return the number of bins along the second axis.
     */
    size_t m( void ) { return( _m ); }

    /*! \brief Return the step size along along the first axis.
     */
    double nstep( void ) { return( _nstep ); }

    /*! \brief Return the step size along along the second axis.
     */
    double mstep( void ) { return( _mstep ); }

    /*! \brief Return the coordinate along the first axis on bin \a i.
     */
    double icoord( size_t i ) { return( _range[0] + i*(_range[2]-_range[0]) / (_n-1.0) ); }

    /*! \brief Return the coordinate along the second axis on bin \a j.
     */
    double jcoord( size_t j ) { return( _range[1] + j*(_range[3]-_range[1]) / (_m-1.0) ); }

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
    void get_range( double range[4] ) { 
	range[0] =  _range[0];
	range[1] =  _range[1];
	range[2] =  _range[2];
	range[3] =  _range[3];
    }
    
    /*! \brief Return bin range.
     *
     *  Returns minimum and maximum values on any bin in histogram.
     */
    void get_bin_range( double &min, double &max );
    
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













