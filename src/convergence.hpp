/*! \file convergence.hpp
 *  \brief Vlasov system convergence tester
 */

/* Copyright (c) 2011 Taneli Kalvas. All rights reserved.
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
 * taneli.kalvas@jyu.fi.
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

#ifndef CONVERGENCE_HPP
#define CONVERGENCE_HPP 1


#include <iostream>
#include <stdint.h>
#include <vector>
#include "scalarfield.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"


/*! \brief Vlasov system convergence tester.
 */
class Convergence {

    int32_t                 _iter;              /*!< \brief Iteration counter. */

    std::vector<double>     _epot_hist;         /*!< \brief Convergence history of epot. */
    double                  _epot_lim;          /*!< \brief Convergence limit. */
    double                  _epot_absf;         /*!< \brief Absolute error coefficient. */
    double                  _epot_relf;         /*!< \brief Relative error coefficient. */
    MeshScalarField        *_epot_old;          /*!< \brief Epot of last round. */
    const MeshScalarField  *_epot;              /*!< \brief Epot pointer. */

    std::vector<double>     _scharge_hist;      /*!< \brief Convergence history of scharge. */
    double                  _scharge_lim;       /*!< \brief Convergence limit. */
    double                  _scharge_absf;      /*!< \brief Absolute error coefficient. */
    double                  _scharge_relf;      /*!< \brief Relative error coefficient. */
    MeshScalarField        *_scharge_old;       /*!< \brief Scharge of last round. */
    const MeshScalarField  *_scharge;           /*!< \brief Scharge pointer. */

    struct EmitPoint {
	double              _epsilon;           /*!< \brief Epsilon. */
	double              _alpha;             /*!< \brief Alpha. */
	double              _beta;              /*!< \brief Beta. */
	double              _xave;              /*!< \brief x average. */
	double              _xpave;             /*!< \brief x' average. */

	EmitPoint();
	EmitPoint( double x );
	EmitPoint( const Emittance &emit );

	double &operator[]( int i );
    };

    std::vector<double>     _emit_hist[2];      /*!< \brief Convergence history of emittance. */
    double                  _emit_lim;          /*!< \brief Convergence limit. */
    double                  _emit_absf;         /*!< \brief Absolute error coefficient. */
    double                  _emit_relf;         /*!< \brief Relative error coefficient. */
    EmitPoint               _emit_old[2];       /*!< \brief Emittances of last round. */
    const ParticleDataBase *_emit_pdb;          /*!< \brief Particle database pointer. */
    coordinate_axis_e       _emit_axis;         /*!< \brief Emittance plane axis. */
    double                  _emit_val;          /*!< \brief Emittance plane location. */
    const Emittance        *_emit[2];           /*!< \brief Emittance pointers. */
   
public:

    /*! \brief Constructor for convergence class.
     */
    Convergence();

    /*! \brief Destructor for convergence class.
     */
    ~Convergence();

    /*! \brief Evaluate convergence of iteration round.
     *
     *  Increase iteration round counter and calculate convergence
     *  estimates for each followed feature. Returns true if
     *  convergence limits are reached. False is returned otherwise.
     *  Convergence estimates are saved to be analyzed. A summary is
     *  printed to cout if verbose printing is enabled.
     */
    bool evaluate_iteration( void );

    /*! \brief Print the history of convergence to stream.
     */
    void print_history( std::ostream &os ) const;

    /*! \brief Add electric potential to be followed.
     */
    void add_epot( const MeshScalarField &epot, double absf, double relf, double lim );

    /*! \brief Add space charge density to be followed.
     */
    void add_scharge( const MeshScalarField &scharge, double absf, double relf, double lim );

    /*! \brief Add trajectory diagnostics to be followed.
     */
    void add_tdiag( const ParticleDataBase &pdb, coordinate_axis_e axis,
		    double val, double absf, double relf, double lim );

    /*! \brief Add trajectory diagnostics to be followed.
     */
    void add_tdiag( const Emittance *emit1, double absf, double relf, double lim ) {
	add_tdiag( emit1, NULL, absf, relf, lim );
    }

    /*! \brief Add trajectory diagnostics to be followed.
     */
    void add_tdiag( const Emittance *emit1, const Emittance *emit2, 
		    double absf, double relf, double lim );

    /*! \brief Set emittance for trajectory diagnostics.
     */
    void set_emittance( const Emittance *emit1, const Emittance *emit2 = NULL );

};


#endif
