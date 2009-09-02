/*! \file problem.hpp
 *  \brief Header file defining abstract base class %Problem.
 */


#ifndef PROBLEM_HPP
#define PROBLEM_HPP 1


#include "crowmatrix.hpp"
#include "mvector.hpp"


/*! \brief Abstract base class for linear/non-linear problem, which
 *  can be described by a non-linear system of equations.
 */
class Problem {

public:

    /*! \brief Virtual destructor.
     */
    virtual ~Problem() {}

    /*! \brief Return const pointers to the matrix \a A and vector \a
     *  B of the linear problem.
     */
    virtual void get_vecmat( const Matrix **A, const Vector **B ) const = 0;

    /*! \brief Return const pointers to jacobian matrix and residual
     *  vector of the problem to \a J and \a R at \a X.
     */
    virtual void get_resjac( const Matrix **J, const Vector **R, const Vector &X ) const = 0;

    /*! \brief Return true if problem is linear.
     */
    virtual bool linear( void ) const = 0;

};


#endif














