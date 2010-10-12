/*! \file epot_problem.hpp
 *  \brief Header file for epot_problem.hpp
 */

/* Copyright (c) 2005-2010 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * tvkalvas@cc.jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */

#ifndef EPOT_PROBLEM_HPP
#define EPOT_PROBLEM_HPP 1


#include <iostream>
#include <stdint.h>
#include "problem.hpp"
#include "solver.hpp"
#include "scalarfield.hpp"
#include "geometry.hpp"
#include "vec3d.hpp"


/*! \brief Plasma modes
 *
 *  Selection of modes for plasma calculation in electric potential
 *  problem. In a typical positive ion plasma extraction simulation
 *  the plasma mode is PLASMA_PEXP_INITIAL in the first iteration
 *  round to use initial guess for plasma meniscus
 *  location. Thereafter PLASMA_PEXP is used. For negative ion
 *  extraction the PLASMA_NSIMP_INITIAL is used for the first iteration
 *  and PLASMA_NSIMP thereafter. PLASMA_INITIAL is a macro, which
 *  equals to PLASMA_PEXP_INITIAL. It exists for backward compatibility.
 *
 */
enum plasma_mode_e {PLASMA_NONE = 0, PLASMA_PEXP_INITIAL, PLASMA_PEXP, 
		    PLASMA_NSIMP_INITIAL, PLASMA_NSIMP};

#define PLASMA_INITIAL PLASMA_PEXP_INITIAL

/*! \brief Class for constructing the linear/nonlinear problem for the
    solver.
 *
 *  %EpotProblem class constructs the Poisson (equation) problem in finite
 *  difference form from Geometry (mesh) and various parameters and it
 *  presents the problem to the Solver via matrix/vector
 *  representation. In case of linear problem (no plasma model), the
 *  Poisson equation 
 *  \f[ \nabla^2 \phi = -\frac{\rho}{\epsilon_0} \f]
 *  is
 *  \f[ \frac{\partial^2 \phi}{\partial x^2} = -\frac{\rho}{\epsilon_0} \f]
 *  in 1D, which is discretized into 
 *  \f[ \phi_{i-1} + \phi_{i+1} - 2\phi_{i} = -h^2 \frac{\rho_{i}}{\epsilon_0} \f]
 *  using finite differences. In 2D coordinates the discretized form is
 *  \f[ \phi_{i-1,j} + \phi_{i+1,j} + \phi_{i,j-1} + \phi_{i,j+1} - 4\phi_{i,j}
    = -h^2 \frac{\rho_{i,j}}{\epsilon_0} \f]
 *  and in 3D it is
 *  \f[ \phi_{i-1,j,k} + \phi_{i+1,j,k} + \phi_{i,j-1,k} + \phi_{i,j+1,k} 
    + \phi_{i,j,k-1} + \phi_{i,j,k+1} - 6\phi_{i,j} = -h^2 \frac{\rho_{i,j,k}}{\epsilon_0}. \f]
 *  In cylindrical coordinates the Poisson equation is
 *  \f[ \frac{\partial^2 \phi}{\partial r^2} + \frac{1}{r} \frac{\partial \phi}{\partial r}
    + \frac{1}{r^2} \frac{\partial^2 \phi}{\partial \theta^2} + \frac{\partial^2 \phi}{\partial z^2} 
    = -\frac{\rho}{\epsilon_0}, \f]
 *  where \f$ \frac{\partial^2 \phi}{\partial \theta^2} = 0 \f$ because of
 *  cylindrical symmetry of the simulations. Therefore the discretized form
 *  becomes
 *  \f[ \phi_{i-1,j} + \phi_{i+1,j} 
    + \left( 1 - \frac{h}{2r_j} \right) \phi_{i,j-1} 
    + \left( 1 + \frac{h}{2r_j} \right) \phi_{i,j+1} - 4\phi_{i,j}
    = -h^2 \frac{\rho_{i,j}}{\epsilon_0}, \f]
 *  where \f$ \frac{h}{2r_j} = \frac{h}{2jh} = \frac{1}{2j} \f$ because the 
 *  radius \f$ r_j = hj \f$. At the symmetry axis there is an exception because both 
 *  \f$ r \f$ and \f$ \frac{\partial \phi}{\partial r} \f$ approach zero. By 
 *  using Bernoulli-L'Hopital rule we can evaluate
 *  \f[ \lim_{r \rightarrow 0} \frac{1}{r} \frac{\partial \phi}{\partial r} 
    = \lim_{r \rightarrow 0} \frac{\partial^2 \phi}{\partial r^2}, \f]
 *  which turns the Poisson equation to
 *  \f[ 2 \frac{\partial^2 \phi}{\partial r^2} + \frac{\partial^2 \phi}{\partial z^2} 
    = -\frac{\rho}{\epsilon_0} \f]
 *  on axis. Discretation of the equation gives us
 *  \f[ \phi_{i-1,j} + \phi_{i+1,j} 
    + 2\phi_{i,j+1} + 2\phi_{i,j-1} - 6\phi_{i,j}
    = -h^2 \frac{\rho_{i,j}}{\epsilon_0}. \f]
 *  Here \f$ \phi_{i,j-1} = \phi_{i,j+1} \f$ so that the final form is
 *  \f[ \phi_{i-1,j} + \phi_{i+1,j} 
    + 4\phi_{i,j+1} - 6\phi_{i,j}
    = -h^2 \frac{\rho_{i,j}}{\epsilon_0}. \f]
 *
 *  In addition to the Poisson equation the problem matrix and vector
 *  also contain finite difference representations of the boundary
 *  conditions.  The Dirichlet boundary condition is defined by
 *  constant potential at boundary, i.e.  \f$ \phi_{i,j} =
 *  \phi_{\mathrm{const}} \f$. The Neumann boundary condition can be
 *  defined as first order discretation 
 *  \f[ \frac{\phi_{i+1}-\phi_{i}}{h} = N_{\mathrm{const}} \f]
 *  or second order discretation
 *  \f[ \frac{-\phi_{i+2}+4\phi_{i+1}-3\phi_{i}}{2h} = N_{\mathrm{const}} \f]
 *  selected by the user.
 *
 *  The plasma problems are described by using nonlinear models for
 *  screening charges in the plasma. For positive ion extraction for
 *  example, the screening charge is an electron population at the
 *  plasma potential \f$ \phi_p \f$ with a thermal energy distribution
 *  with temperature \f$ T_e \f$. The screening charge is therefore
 *  \f[ \rho_e = \rho_{e0} \exp \left( \frac{\phi-\phi_p}{kT_e/e} \right), \f]
 *  where electron charge density at plasma potential \f$ (\rho_{e0}) \f$ 
 *  is the same as the total positive beam space charge density for enabling
 *  plasma neutrality.
 */
class EpotProblem : public Problem {

    /*! \brief Class nodes to degrees of freedom mapping.
     *
     *  Uses running numbers starting from 0 to point to vectors with
     *  free variables (electric potential and matrices during solving
     *  for example). If a node is solid interior point or dirichlet
     *  boundary, a fixed negative number corresponding to the solid
     *  is used, i.e. from -1 to -6 for boundaries and starting from
     *  -7 for electrodes.
     */
    class Node2DoF {
	Int3D         _size;          /*!< \brief Size of mesh */
	int32_t      *_n2d;           /*!< \brief Nodes to degrees of freedom array. */
	
    public:
	
	Node2DoF() : _size(0), _n2d(0) {}
	Node2DoF( Int3D size ) : _size(size) {
	    _n2d = new int32_t[_size[0]*_size[1]*_size[2]];
	}
	~Node2DoF() { delete _n2d; }
	
	void resize( Int3D size ) {
	    _size = size;
	    if( _n2d )
		delete _n2d;
	    _n2d = new int32_t[_size[0]*_size[1]*_size[2]];
	}
	
	int32_t &operator()( int i ) 
	    { return( _n2d[i] ); }
	int32_t &operator()( int i, int j ) 
	    { return( _n2d[i+j*_size[0]] ); }
	int32_t &operator()( int i, int j, int k ) 
	    { return( _n2d[i+j*_size[0]+k*_size[0]*_size[1]] ); }
	
	const int32_t &operator()( int i ) const
	    { return( _n2d[i] ); }
	const int32_t &operator()( int i, int j )  const 
	    { return( _n2d[i+j*_size[0]] ); }
	const int32_t &operator()( int i, int j, int k ) const 
	    { return( _n2d[i+j*_size[0]+k*_size[0]*_size[1]] ); }
	
	void debug_print( void ) const; 
    };

    int32_t             _nodecount;     /*!< \brief Number of nodes. */
    int32_t             _dof;           /*!< \brief Degrees of freedom. */
    Node2DoF            _n2d;           /*!< \brief Nodes to degrees of freedom map. */
    CRowMatrix         *_fd_mat;        /*!< \brief Finite Difference (linear) matrix. */
    Vector             *_fd_vec;        /*!< \brief Finite Difference (linear) vector. */

    const Geometry     *_g;             /*!< \brief Geometry for solve in process. */
    mutable CRowMatrix *_fd_mat2;       /*!< \brief Working matrix for solver. */
    mutable Vector     *_fd_vec2;       /*!< \brief Working vector 1 for solver. */
    mutable Vector     *_fd_vec3;       /*!< \brief Working vector 2 for solver. */

    int32_t             _neumann_order; /*!< \brief Neumann boundary order (1 or 2). */
    bool                _smooth_solid;  /*!< \brief Enable smooth solids. */

    plasma_mode_e       _plasma;        /*!< \brief Plasma simulation mode. */

    double              _rhoe;          /*!< \brief Electron charge density (C/m3), < 0. */
    double              _Te;            /*!< \brief Electron thermal energy, > 0. */
    double              _Up;            /*!< \brief Plasma potential, > 0. */

    std::vector<double> _rhoi;          /*!< \brief Charge density for positive ions, 
					 *   first fast protons, then thermal ions */
    std::vector<double> _Ei;            /*!< \brief Energy for positive ions,
					 *   first fast protons, then thermal ions */

    double              _meniscus_x;    /*!< \brief Initial plasma meniscus x-coordinate. */
    int32_t             _meniscus_i;    /*!< \brief Initial plasma meniscus i-coordinate (calculated). */

    Solver             *_solver;        /*!< \brief Solver for solving problem.. */


    void set_link( CRowMatrix &A, Vector &B, 
		   int32_t a, int32_t b, double val );

    void add_initial_plasma( int32_t i, int32_t j, int32_t k, 
			     CRowMatrix &A, Vector &B, Node2DoF &n2d );

    void add_vacuum_node( int32_t i, int32_t j, int32_t k, 
			  CRowMatrix &A, Vector &B, Node2DoF &n2d );

    void add_neumann_node( signed char a, int32_t i, int32_t j, int32_t k, 
			   CRowMatrix &A, Vector &B, Node2DoF &n2d );

    void add_solid_edge_node( signed char a, int32_t i, int32_t j, int32_t k, 
			      CRowMatrix &A, Vector &B, Node2DoF &n2d );
    
    void clear_problem( void );

public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Default constructor.
     */
    EpotProblem();

    /*! \brief Constructor for loading problem from a file.
     */
    EpotProblem( std::istream &s );

    /*! \brief Destructor for problem.
     */
    ~EpotProblem();
    
/* ************************************** *
 * Problem constructing and solving       *
 * ************************************** */

    /*! \brief Set Neumann boundary order.
     *
     *  Valid values are 1 and 2 (default).
     */
    void set_neumann_order( int32_t order );

    /*! \brief Enable smooth solid edges
     *
     *  Smooth edges are enabled by default.
     */
    void enable_smooth_solids( bool enable );

    /*! \brief Define initial plasma to the problem.
     */
    void set_initial_plasma( double Up, double meniscus_x );

    /*! \brief Enable plasma model for positive ion extraction problem.
     */
    void set_pexp_plasma( double rhoe, double Te, double Up );

    /*! \brief Define initial plasma boundary location to negative ion
     *  extraction problem.
     */
    void set_nsimp_initial_plasma( double meniscus_x );

    /*! \brief Enable plasma model for negative ion extraction problem.
     *
     *  The positive (analytic) space charges for the negative ion
     *  plasma extraction are set using this function. The positive
     *  ions consist of fast (directed) protons and any number of
     *  thermal positive ions trapped at the plasma boundary in the
     *  zero potential.
     *
     *  The parameters set are \a rhop, the space charge density of
     *  protons and \a Ep, the energy of protons at zero
     *  potential. Vectors \a rhoi and \a Ei are used to set the space
     *  charge densities and thermal energies of the trapped ions.
     */
    void set_nsimp_plasma( double rhop, double Ep, 
			   std::vector<double> rhoi, std::vector<double> Ei );

    /*! \brief Construct matrix form of the problem.
     *
     *  Requires that mesh is build for geometry \a g.
     */
    void construct( const Geometry &g );

    /*! \brief Set solver to be used for the problem.
     */
    void set_solver( Solver &s );

    /*! \brief Solve the problem.
     *
     *  The \a epot field is used as an initial guess for the
     *  solver. The space charge density field \a scharge is added to
     *  the problem vector before solving. The solution is returned in
     *  \a epot.
     */
    void solve( ScalarField &epot, const ScalarField &scharge ) const;

/* ************************************** *
 * Solver interface                       *
 * ************************************** */

    /*! \brief Return const pointers to the matrix \a A and vector \a
     *  B of the linear problem.
     *
     *  This function should only be called by Solver, while method
     *  Problem::solve() is running and only in case of linear
     *  problem.
     */
    void get_vecmat( const Matrix **A, const Vector **B ) const;

    /*! \brief Return const pointers to jacobian matrix and residual
     *  vector of the problem to \a J and \a R at \a X.
     *
     *  This function should only be called by Solver, while method
     *  Problem::solve() is running and only in case of nonlinear
     *  problem.
     */
    void get_resjac( const Matrix **J, const Vector **R, const Vector &X ) const;

    /*! \brief Return true if problem is linear.
     */
    bool linear( void ) const;

/* ************************************** *
 * Misc                                   *
 * ************************************** */

    /*! \brief Get degrees of freedom.
     */
    int get_dof( void ) const { return( _dof ); }

    /*! \brief Prints internal data to std::cout.
     */
    void debug_print( void ) const;

    /*! \brief Saves problem data to stream.
     */
    void save( std::ostream &s ) const;
};


#endif
















