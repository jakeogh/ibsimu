/*! \file polysolver.hpp
 *  \brief Header file for polynomial solvers
 */


#ifndef POLYSOLVER_HPP
#define POLYSOLVER_HPP 1


#include <stdint.h>


/*! \brief Solve quadric equation a*x^2 + b*x + c = 0.
 *
 *  Solves the quadric equation. Can also handle linear case if a is
 *  zero. Returns the number of roots found. Roots are filled in
 *  ascending order.
 */
uint32_t solve_quadratic( double a, double b, double c, 
			  double *x0, double *x1 );


/*! \brief Solve cubic equation a*x^3 + b*x^2 + c*x + d = 0.
 *
 *  Solves the cubic equation.  Can also handle lower order
 *  polynomials if coefficients are zero. Returns the number of
 *  roots found. Roots are filled in ascending order.
 */
uint32_t solve_cubic( double a, double b, double c, double d, 
		      double *x0, double *x1, double *x2 );


/*! \brief Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d = 0.
 *
 *  Solves the quartic equation. Returns the number of roots
 *  found. Roots are filled in ascending order.
 */
uint32_t solve_quartic( double a, double b, double c, double d,
			double *x0, double *x1, double *x2, double *x3 );


#endif













