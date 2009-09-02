/*! \file qrandom.hpp
 *  \brief Header file defining quasi random number generator class.
 */


#ifndef QRANDOM_HPP
#define QRANDOM_HPP 1


#include <gsl/gsl_qrng.h>
#include <gsl/gsl_roots.h>


/*! \brief Quasi random number generator for N dimensions.
 *
 *  This QRNG can produce quasi random numbers in N independent
 *  dimensions. QRNG includes functions to return uniformly
 *  distributed numbers between 0 and 1 and numbers from a gaussian
 *  distribution.
 */
class QRandom {
    int                 N;        /*!< \brief Number of dimensions for QRNG. */
    gsl_qrng           *qrng;     /*!< \brief Random number generator from gsl. */
    gsl_function_fdf    fdf;      /*!< \brief Function to solve for gaussian transformation. */
    gsl_root_fdfsolver *solver;   /*!< \brief Solver for gaussian transformation. */

    QRandom( const QRandom &qrng ) {}

    /*! \brief Returns number from gaussian distribution, transformed
     *  from uniformly distributed R, where 0 <= R <= 1.
     */
    double transform_gaussian( double R );

public:

    /*! \brief Constructor for QRNG in \a N independent dimensions.
     */
    QRandom( int n );

    /*! \brief Destructor.
     */
    ~QRandom();

    /*! \brief Get gaussian random numbers.
     *
     *  Get next sampling from gaussian quasi random number generator
     *  to \a x with standard deviation of 1. Array \a x must have
     *  space for \a N numbers. */
    void get_gaussian( double *x );
    
    /*! \brief Get uniform and gaussian random numbers.
     * 
     *  Get next sampling from quasi random number generator to \a x
     *  with dimensions marked (to true) in list gaussian mapped
     *  to gaussian distribution with standard deviation of 1. Rest of
     *  the dimensions are linear. Array x must have space for \a N
     *  numbers. */
    void get_part_gaussian( bool *gaussian, double *x );

    /*! \brief Get uniform random numbers.
     * 
     *  Get next sampling from quasi random number generator to \a x.
     *  Array x must have space for \a N numbers. */
    void get( double *x );
};


#endif













