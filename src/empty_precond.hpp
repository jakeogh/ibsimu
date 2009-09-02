/*! \file empty_precond.hpp
 *  \brief Header file for the empty preconditioner class.
 */


#ifndef EMPTY_PRECOND_HPP
#define EMPTY_PRECOND_HPP 1


#include "precond.hpp"


/*! \brief Empty preconditioner class.
 */
class Empty_Precond : public Precond {
public:

    /*! \brief Constructor.
     */
    Empty_Precond() {}

    /*! \brief Destructor.
     */
    ~Empty_Precond() {}

    /*! \brief Return \a x = \a b.
     */
    inline void solve( Vector &x, const Vector &b ) const {
	x = b;
    }
};


#endif














