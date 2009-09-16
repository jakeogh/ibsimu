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
 *  solver library doing LU decomposition to solve linear matrix
 *  problems. This solver includes both linear and nonlinear
 *  solvers. Nonlinear solver based on Newton-Raphson iteration.
 */
class UMFPACKSolver : public Solver {

    void    *_numeric;      /*!< \brief Numeric data for LU decomposition. */

    double   _newton_Reps;  /*!< \brief Accuracy request for Newton-Raphson residual. */
    double   _newton_dXeps; /*!< \brief Accuracy request for Newton-Raphson step. */
    double   _newton_imax;  /*!< \brief Maximum number of Newton-Raphson iterations. */

    static void umfpack_error( const std::string func, int status );

    void umfpack_decompose( const CColMatrix &mat );
    void umfpack_solve( const CColMatrix &mat, const Vector &rhs, Vector &sol,
			bool force_decomposition = false );

public:

    /*! \brief Constructor.
     */
    UMFPACKSolver( double newton_Reps = 1.0e-5, double newton_dXeps = 1.0e-6, 
		   uint32_t newton_imax = 10 );

    /*! \brief Destructor.
     */
    ~UMFPACKSolver();

    /*! \brief Solve problem \a p defined in geometry \a g.  Initial
     *  guess and solution are in vector \a X.
     */
    virtual void solve( const Problem &p, Vector &X );

    /*! \brief Reset solver.
     *
     *  This is a signal from the problem that the problem has changed
     *  and internal caches (if they exist) in the solver should be
     *  resetted.
     */
    virtual void reset( void );

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
};


#endif
