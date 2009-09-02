/*! \file ilu0_precond.hpp
 *  \brief Header file for the zero fill-in incomplete LU preconditioner class.
 */


#ifndef ILU0_PRECOND_HPP
#define ILU0_PRECOND_HPP 1


#include "matrix.hpp"
#include "crowmatrix.hpp"
#include "ccolmatrix.hpp"
#include "precond.hpp"


/*! \brief Zero fill-in incomplete LU preconditioner class.
 */
class ILU0_Precond : public Precond {
    Matrix *_L, *_U;

public:

    /*! \brief Constructor for an ILU0 preconditioner for matrix \a A.
     */
    ILU0_Precond( const Matrix &A );

    /*! \brief Destructor.
     */
    ~ILU0_Precond();

    /*! \brief Prints the values of all internal data to std::cout.
     */
    void debug_print( void ) const;

    /*! \brief Returns a pointer to the internal L matrix.
     */
    const Matrix *get_L( void ) const { return( _L ); }
    
    /*! \brief Returns a pointer to the internal U matrix.
     */
    const Matrix *get_U( void ) const { return( _U ); }

    /*! \brief Solve \a M* \a x = \a b and return \a x.
     */
    void solve( Vector &x, const Vector &b ) const;
};




#endif














