/*! \file trajectory.hpp
 *  \brief Header file defining particle trajectory helpers
 */


#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP 1


enum trajectory_rep_e {
    TRAJ_EMPTY = 0,
    TRAJ_LINEAR,
    TRAJ_QUADRATIC,
    TRAJ_CUBIC
};


/*! \brief Trajectory representation between two calculated points in 1d.
 *
 *  Uses either linear-linear, quadratic-linear or cubic-quadratic
 *  representation for location \a x and velocity \a v. Time is
 *  presented as a parametric, scaled variable \a K ranging from 0 to 1.
 */
class TrajectoryRep1D {

    trajectory_rep_e _rep;
    double           _A, _B, _C, _D, _E;

public:

    /*! \brief Default constructor for empty representation.
     */
    TrajectoryRep1D() : _rep(TRAJ_EMPTY) {}

    /*! \brief Constructor for representation of trajectory from \a
     *  (x1,v1) to \a (x2,v2) in time \a dt.
     *
     *  Can be forced to use a specified representation type by
     *  setting \a force. Defaults to TRAJ_EMPTY, which means that the
     *  highest numerically stable polynomial is automatically used.
     */
    TrajectoryRep1D( double dt, double x1, double v1, double x2, double v2,
		     trajectory_rep_e force = TRAJ_EMPTY );

    /*! \brief Destructor.
     */
    ~TrajectoryRep1D() {}

    /*! \brief Construct representation of trajectory from \a (x1,v1)
     *  to \a (x2,v2) in time \a dt.
     *
     *  Can be forced to use a specified representation type by
     *  setting \a force. Defaults to TRAJ_EMPTY, which means that the
     *  highest numerically stable polynomial is automatically used.
     */
    void construct( double dt, double x1, double v1, double x2, double v2,
		    trajectory_rep_e force = TRAJ_EMPTY );

    /*! \brief Calculate location \a x and velocity \a v at parametric
     *  time \a K.
     */
    void coord( double &x, double &v, double K );

    /*! \brief Solves for trajectory intersection with location.
     *
     *  Solves the trajectory intersection with location \a x. Saves
     *  the valid solutions to array \a K in increasing order and
     *  returns the number of solutions saved. The allowed limit for
     *  parametric time \a K is 0 < K =< 1, if \a include_start is
     *  false or 0 <= K =< 1, if \a include_start is true.
     */
    int solve( double K[3], double x );

    /*! \brief Prints internal data to std::cout.
     */
    void debug_print( void ) const;
};


#endif













