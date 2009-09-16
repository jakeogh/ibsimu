/*! \file scalarfield.hpp
 *  \brief Header file defining scalar field class.
 */


#ifndef SCALARFIELD_HPP
#define SCALARFIELD_HPP 1


#include <iostream>
#include "geometry.hpp"
#include "vec3d.hpp"


/*! \page p_epot %Electric potential field
 *
 *  The electric potential field is defined using a ScalarField class
 *  object, the same which is used for \subpage p_scharge
 *  "space charge". The electric potential field contains the potential
 *  values in a rectangular mesh. In a regular simulation the mesh is
 *  defined by the Geometry class. The electric potential field can be
 *  calculated in IBSimu using the mesh, solid and boundary condition
 *  information from a %Geometry object by defining a electric
 *  potential problem EpotProblem. This problem can the be solved with
 *  suitable solvers in the Ion Beam Simulation Package.
 *
 *  The electric potential field can be used to define an electric
 *  field using EpotEfield class, which provides electric field values
 *  by differentiating and interpolating the potential field
 *  on-line. The whole chain from the geometry to a solved electric
 *  field is shown by the following (partial) code example:
 *
\code
#include <geometry.hpp>
#include <func_solid.hpp>
#include <scalarfield.hpp>
#include <epot_problem.hpp>
#include <epot_efield.hpp>
#include <umfpack_solver.hpp>

int main( void )
{
    ...

    EpotProblem problem;
    problem.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );

    UMFPACKSolver solver;
    problem.set_solver( solver );
    problem.solve( epot, scharge );

    EpotEfield ef( geom, epot );

    return( 0 );
}
\endcode
 *
 */


/*! \brief Scalar field class.
 *
 *  %ScalarField is intended to be used as a container for scalar
 *  field data. The IBSimu Package uses scalar fields for storing
 *  electric potential and space charge density fields. A separate
 *  Efield exists for electric fields. The scalar field class
 *  provides a subset of vector operations to copy, sum and scale
 *  scalar fields. The class also includes operators for indexed
 *  access to mesh elements and interpolation functions for linearly
 *  interpolated smooth field data. The interpolation function uses
 *  linear extrapolation of field outside the defined mesh. This way
 *  it can be ensured that %ScalarField returns sensible values even
 *  close by to the edges of the geometry.
 */
class ScalarField {
    geom_mode_e _geom_mode; /*!< \brief Geometry mode */
    Int3D       _size;      /*!< \brief Size of mesh */
    Vec3D       _origo;     /*!< \brief Location of mesh point (0,0,0) [m] */
    Vec3D       _max;       /*!< \brief Location of mesh point (size[0]-1,size[1]-1,size[2]-1) [m] */
    double      _h;         /*!< \brief Length of mesh step [m] */
    double      _div_h;     /*!< \brief One over length of mesh step [1/m] */
    double     *_F;         /*!< \brief Scalar field data */

public:


    /*! \brief Default constructor.
     */
    ScalarField() : _geom_mode(MODE_3D), _h(1.0), _div_h(1.0), _F(0) {}

    /*! \brief Constructor for geometry from \a g.
     *
     *  Returns a new scalar field with geometry parameters (including
     *  mesh size) set from \a g. The field is set to zero in all
     *  locations.
     */
    ScalarField( const Geometry &g );

    /*! \brief Constructor for set geometry.
     *
     *  Returns a new scalar field with geometry set according to
     *  parameters. The field is set to zero in all locations.
     */
    ScalarField( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h );

    /*! \brief Constructor for loading scalar field from a file.
     */
    ScalarField( std::istream &s );

    /*! \brief Copy constructor.
     */
    ScalarField( const ScalarField &f );

    /*! \brief Destructor.
     */
    ~ScalarField();

    /*! \brief Returns geometry mode.
     */
    inline geom_mode_e geom_mode( void ) const { return( _geom_mode ); }

    /*! \brief Returns size array of geometry.
     */
    inline Int3D size( void ) const { return( _size ); }

    /*! \brief Returns size of solid mesh in direction \a i.
     */
    inline int32_t size( int i ) const { return( _size[i] ); }
   
    /*! \brief Returns number of nodes in the mesh.
     */
    inline int32_t nodecount( void ) const { return( _size[0]*_size[1]*_size[2] ); }

    /*! \brief Returns origo vector of geometry.
     */
    inline Vec3D origo( void ) const { return( _origo ); }

    /*! \brief Returns \a i-th component of vector origo.
     */
    inline double origo( int i ) const { return( _origo[i] ); }

    /*! \brief Returns vector pointing to the last mesh point opposite
     *  of origo.
     */
    Vec3D max( void ) const { return( _max ); }

    /*! \brief Returns \a i-th component of vector pointing to the
     *  last mesh point opposite of origo.
     */
    double max( int i ) const { return( _max[i] ); }

    /*! \brief Returns mesh cell size.
     */
    inline double h( void ) const { return( _h ); }

    /*! \brief Returns inverse mesh cell size.
     */
    inline double div_h( void ) const { return( _div_h ); }

    /*! \brief Clears the field.
     */
    void clear();

    /*! \brief Resets the field geometry.
     *
     *  Sets the field geometry according to the parameters and clears
     *  the field to zero in all locations.
     */
    void reset( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h );

    /*! \brief Search minimum and maximum values of scalar field.
     *
     *  It should be noted that %ScalarField object itself doesn't
     *  know about virtual potentials at electrode edges and therefore
     *  this function doesn't know how to separate these points from
     *  other points. For electric potential fields use
     *  epot_get_minmax() instead.
     */
    void get_minmax( double &min, double &max ) const;

    /*! \brief Search minimum and maximum electric potential values of
     *  electric potential field.
     *
     *  Searches the minimum and maximum electric potential field
     *  values skipping the virtual potentials of electrode edges
     *  defined by the geometry and adding in the electrode potentials.
     */
    void epot_get_minmax( const Geometry &g, double &min, double &max ) const;

    /*! \brief Copy operator.
     */
    ScalarField &operator=( const ScalarField &f );

    /*! \brief Accumulation operator.
     */
    ScalarField &operator+=( const ScalarField &f );

    /*! \brief Scaling operator.
     */
    ScalarField &operator*=( double x );

    /*! \brief Inverse scaling operator.
     */
    ScalarField &operator/=( double x );

    /*! \brief Operator for pointing to elements of the field mesh.
     */
    inline double &operator()( int32_t i ) { 
	return( _F[i] ); }

    /*! \brief Operator for pointing to elements of the field mesh.
     */
    inline double &operator()( int32_t i, int32_t j ) {
	return( _F[i + j*_size[0]] ); }

    /*! \brief Operator for pointing to elements of the field mesh.
     */
    inline double &operator()( int32_t i, int32_t j, int32_t k ) {
	return( _F[i + (j + k*_size[1])*_size(0)] ); }

    /*! \brief Operator for pointing to elements of the field mesh.
     */
    inline const double &operator()( int32_t i ) const { 
	return( _F[i] ); }

    /*! \brief Operator for pointing to elements of the field mesh.
     */
    inline const double &operator()( int32_t i, int32_t j ) const {
	return( _F[i + j*_size[0]] ); }

    /*! \brief Operator for pointing to elements of the field mesh.
     */
    inline const double &operator()( int32_t i, int32_t j, int32_t k ) const {
	return( _F[i + (j + k*_size[1])*_size(0)] ); }

    /*! \brief Operator for getting linearly interpolated field values.
     *
     *  The field is interpolated linearly to get the field value at
     *  \a x. If \a x is outside the mesh, the field is extrapolated
     *  linearly using the field points to \a x. This provides correct
     *  field values also close to the mesh boundaries.
     */
    double operator()( Vec3D x ) const;

    /*! \brief Saves scalar field data to stream.
     */
    void save( std::ostream &s ) const;

    /*! \brief Prints internal data to std::cout.
     */
    void debug_print( void ) const;
};


#endif













