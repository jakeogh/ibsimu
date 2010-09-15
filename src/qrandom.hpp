/*! \file qrandom.hpp
 *  \brief Header file for qrandom.hpp
 */

/* Copyright (c) 2005-2009 Taneli Kalvas. All rights reserved.
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
 * tvkalvas@cc.jyu.fi.
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















