/*! \file bicgstab.hpp
 *  \brief Header file for BiCGSTAB solver.
 */


#ifndef BICGSTAB_HPP
#define BICGSTAB_HPP 1


#include "matrix.hpp"
#include "precond.hpp"


/*! \brief BiCGSTAB linear matrix equation solver.
 *
 *  Solves linear matrix equations of form A*X=B.
 *
 *  \return True on success, false on errors
 *  \param mat Matrix on the linear equation
 *  \param rhs Right hand side vector on the linear equation
 *  \param sol Solution vector on the linear equation
 *  \param pc Preconditioner used by the iterator
 *  \param imax Maximum number of iterations to be done in input, 
 *   number of iterations actually done on output
 *  \param eps Accuracy request for iteration on input, 
 *   estimated error achieved on output
 */
bool bicgstab( const Matrix &mat, const Vector &rhs, Vector &sol,
	       const Precond &pc, int &imax, double &eps );



#endif
















