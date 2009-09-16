/*! \file gs_solver.hpp
 *  \brief Header file defining class GSSolver.
 */


#ifndef GS_SOLVER_HPP
#define GS_SOLVER_HPP 1


#include <iostream>
#include "solver.hpp"
#include "geometry.hpp"
#include "problem.hpp"
#include "matrix.hpp"
#include "crowmatrix.hpp"
#include "mvector.hpp"


/*! \brief Gauss-Seidel based solver implementation.
 *
 *  An implementation of virtual class Solver using
 *  Gauss-Seidel/Successive over relaxation method for solving matrix
 *  equation. Linear and nonlinear solvers. Nonlinear solver based on
 *  Newton-Raphson iteration.
 *
 */
class GSSolver : public Solver {
    double   _eps;
    int      _imax;
    double   _w;

    double   _newton_Reps;  /*!< \brief Accuracy request for Newton-Raphson residual. */
    double   _newton_dXeps; /*!< \brief Accuracy request for Newton-Raphson step. */
    double   _newton_imax;  /*!< \brief Maximum number of Newton-Raphson iterations. */

    static void gauss_seidel_error( const std::string func, int status );

public:

    /*! \brief Constructor.
     */
    GSSolver( double w = 1.66,
	      double eps = 1.0e-6, 
	      int imax = 10000, 
	      double newton_Reps = 1.0e-5, 
	      double newton_dXeps = 1.0e-6, 
	      int newton_imax = 10 );

    /*! \brief Destructor.
     */
    ~GSSolver() {}

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

    /*! \brief Direct interface to gauss_seidel matrix solver.
     *
     *  Solves matrix problem \a mat * \a sol = \a rhs. The matrix is
     *  sorted to be on the ascending order as required by the GAUSS_SEIDEL
     *  library.
     */
    static bool gauss_seidel( const CRowMatrix &mat, const Vector &rhs, Vector &sol,
			      int &imax, double &eps, double w );
};


#endif
