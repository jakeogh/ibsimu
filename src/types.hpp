/*! \file types.hpp
 *  \brief Header defining generally used enumeration types and corresponding printable strings.
 */


#ifndef TYPES_HPP
#define TYPES_HPP 1


/*! \brief Geometry mode enum.
 *
 *  Simulation geometry can be either 1D (\a MODE_1D), planar 2D (\a
 *  MODE_2D), planar 3D (\a MODE_3D) or it can be defined as cylindrical
 *  symmetrical 2D with coordinates \a x and \a r (\a MODE_CYL).
 */
enum geom_mode_e {
    MODE_1D = 0,
    MODE_2D,
    MODE_CYL, 
    MODE_3D
};


/*! \brief Boundary type enum.
 *
 *  %Solid boundary may have either essential, Dirichlet boundary
 *  condition \a BOUND_DIRICHLET or natural, Neumann boundary
 *  condition \a BOUND_NEUMANN.
 */
enum bound_e {
    BOUND_DIRICHLET = 0,
    BOUND_NEUMANN
};


/*! \brief Coordinate axis identifier.
 */
enum coordinate_axis_e {
    AXIS_X = 0,
    AXIS_Y,
    AXIS_R,
    AXIS_Z
};


extern const char *coordinate_axis_string[];
extern const char *coordinate_axis_string_with_unit[];


/*! \brief Type of diagnostic for trajectories.
 */
enum trajectory_diagnostic_e {
    DIAG_T = 0,
    DIAG_X,
    DIAG_VX,
    DIAG_Y,
    DIAG_R,
    DIAG_VY,
    DIAG_VR,
    DIAG_W,
    DIAG_VTHETA,
    DIAG_Z,
    DIAG_VZ,
    DIAG_CURR,
    DIAG_XP,
    DIAG_YP,
    DIAG_RP,
    DIAG_ZP,
    DIAG_EK,
    DIAG_QM
};


extern const char *trajectory_diagnostic_string[];
extern const char *trajectory_diagnostic_string_with_unit[];






#endif













