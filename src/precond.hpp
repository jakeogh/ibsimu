/*! \file precond.hpp
 *  \brief Header file for the preconditioner class.
 */


#ifndef PRECOND_HPP
#define PRECOND_HPP 1


#include "mvector.hpp"


/*! \brief Abstract base preconditioner class.
 */
class Precond {
public:

    /*! \brief Virtual destructor.
     */
    virtual ~Precond() {};


    /*! \brief Solve \a M* \a x = \a b and return x. Here \a M is a
     *  preconditioner matrix built for matrix \a A. See
     *  implementations of Precond for more information.
     */
    virtual void solve( Vector &x, const Vector &b ) const = 0;
};


#endif














