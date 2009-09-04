/*! \file bicgstab_solver.hpp
 *  \brief Header file defining class BiCGSTABSolver.
 */


#ifndef BICGSTAB_SOLVER_HPP
#define BICGSTAB_SOLVER_HPP 1


#include <iostream>
#include "solver.hpp"
#include "geometry.hpp"
#include "problem.hpp"
#include "matrix.hpp"
#include "mvector.hpp"


/*! \brief BiCGSTAB based solver implementation.
 *
 *  An implementation of virtual class Solver using ILU0
 *  preconditioner and bicgstab() solver function. Linear and
 *  nonlinear solvers.
 *
 */
class BiCGSTABSolver : public Solver {
    double   _eps;          /*!< \brief Accuracy request. */
    int _imax;         /*!< \brief Maximum iteration count. */

    double   _newton_Reps;  /*!< \brief Accuracy request for Newton-Raphson residual. */
    double   _newton_dXeps; /*!< \brief Accuracy request for Newton-Raphson step. */
    double   _newton_imax;  /*!< \brief Maximum number of Newton-Raphson iterations. */

public:

    /*! \brief Constructor.
     */
    BiCGSTABSolver( double eps = 1.0e-6, int imax = 10000,
		    double newton_Reps = 1.0e-5, double newton_dXeps = 1.0e-6, int newton_imax = 10 );

    /*! \brief Destructor.
     */
    ~BiCGSTABSolver() {}

    /*! \brief Solve problem \a p defined in geometry \a g.  Initial
     *  guess and solution are in vector \a X.
     */
    void solve( const Problem &p, Vector &X ) const;

    /*! \brief Sets the accuracy request for BiCGSTAB solver.
     */
    void set_eps( double eps ) {
	_eps = eps;
    }

    /*! \brief Sets maximum iteration count for BiCGSTAB solver.
     */
    void set_imax( int imax ) {
	_imax = imax;
    }

    /*! \brief Sets maximum iteration count for Newton-Raphson steps.
     */
    void set_newton_imax( int newton_imax ) {
	_newton_imax = newton_imax;
    }

    /*! \brief Sets the accuracy request for Newton-Raphson residual.
     */
    void set_newton_residual_eps( int newton_Reps ) {
	_newton_Reps = newton_Reps;
    }

    /*! \brief Sets the accuracy request for Newton-Raphson step size.
     */
    void set_newton_step_eps( int newton_dXeps ) {
	_newton_dXeps = newton_dXeps;
    }

};


#endif















