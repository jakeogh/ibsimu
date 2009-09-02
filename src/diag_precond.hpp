/*! \file diag_precond.hpp
 *  \brief Header file for the diagonal preconditioner class.
 */


#ifndef DIAG_PRECOND_HPP
#define DIAG_PRECOND_HPP 1


#include "matrix.hpp"
#include "precond.hpp"


/*! \brief Diagonal preconditioner class.
 */
class Diag_Precond : public Precond {
    Vector diag;

public:

    /*! \brief Constructor for a diagonal preconditioner for matrix \a A.
     */
    Diag_Precond( const Matrix &A );

    /*! \brief Destructor.
     */
    ~Diag_Precond() {};

    /*! \brief Solve \a M* \a x = \a b and return \a x.
     */
    void solve( Vector &x, const Vector &b ) const;
};




#endif














