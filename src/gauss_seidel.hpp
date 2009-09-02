/*
 *
 *
 */


#ifndef GAUSS_SEIDEL_HPP
#define GAUSS_SEIDEL_HPP 1


#include "crowmatrix.hpp"
#include "mvector.hpp"


bool gauss_seidel( const CRowMatrix &mat, const Vector &rhs, Vector &sol,
		   uint32_t &imax, double &eps, double w );



#endif














