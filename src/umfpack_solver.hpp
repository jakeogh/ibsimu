/*! \file umfpack_solver.hpp
 *  \brief Header file defining class UMFPACKSolver.
 */


#ifndef UMFPACK_SOLVER_HPP
#define UMFPACK_SOLVER_HPP 1


#include <iostream>
#include "solver.hpp"
#include "geometry.hpp"
#include "problem.hpp"
#include "matrix.hpp"
#include "ccolmatrix.hpp"
#include "mvector.hpp"


/*! \brief UMFPACK based solver implementation.
 *
 *  An implementation of virtual class Solver using UMFPACK matrix
 *  solver library. Linear and nonlinear solvers. Nonlinear solver
 *  based on Newton-Raphson iteration.
 *
 */
class UMFPACKSolver : public Solver {
    double   _newton_Reps;  /*!< \brief Accuracy request for Newton-Raphson residual. */
    double   _newton_dXeps; /*!< \brief Accuracy request for Newton-Raphson step. */
    double   _newton_imax;  /*!< \brief Maximum number of Newton-Raphson iterations. */

public:

    /*! \brief Constructor.
     */
    UMFPACKSolver( double newton_Reps = 1.0e-5, double newton_dXeps = 1.0e-6, 
		   uint32_t newton_imax = 10 );

    /*! \brief Destructor.
     */
    ~UMFPACKSolver() {}

    /*! \brief Solve problem \a p defined in geometry \a g.  Initial
     *  guess and solution are in vector \a X.
     */
    void solve( const Problem &p, Vector &X ) const;

    /*! \brief Sets maximum iteration count for Newton-Raphson steps.
     */
    void set_newton_imax( uint32_t newton_imax ) {
	_newton_imax = newton_imax;
    }

    /*! \brief Sets the accuracy request for Newton-Raphson residual.
     */
    void set_newton_residual_eps( uint32_t newton_Reps ) {
	_newton_Reps = newton_Reps;
    }

    /*! \brief Sets the accuracy request for Newton-Raphson step size.
     */
    void set_newton_step_eps( uint32_t newton_dXeps ) {
	_newton_dXeps = newton_dXeps;
    }

    /*! \brief Direct interface to umfpack matrix solver.
     *
     *  Solves matrix problem \a mat * \a sol = \a rhs. The matrix is
     *  sorted to be on the ascending order as required by the UMFPACK
     *  library.
     */
    static void umfpack_solve( CColMatrix &mat, const Vector &rhs, Vector &sol );
};


#endif
