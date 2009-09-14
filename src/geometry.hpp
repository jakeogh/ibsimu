/*! \file geometry.hpp
 *  \brief Header file geometry defining class.
 */


#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP 1


#include <stdint.h>
#include <vector>
#include <iostream>
#include "file.hpp"
#include "vec3d.hpp"
#include "solid.hpp"
#include "types.hpp"


/*! \page p_geometry %Geometry
 *
 *  The simulation geometry is defined as a %Geometry class
 *  object. %Geometry is the most fundamental piece of every
 *  simulation. There cannot be a simulation without a geometry to
 *  define the volume of space used for the simulation. The geometry
 *  class defines the geometry type (either 1D, 2D, cylindrical
 *  symmetry or 3D), a rectangular box in space and an even
 *  rectangular mesh of points inside it for discretizing the
 *  space. This same discretation defined by the %Geometry is used for
 *  electric potential field and space charge field.
 *
 *  The %Geometry class defines the solids or electrodes used in the
 *  simulation. The solids are numbered with running numbers starting
 *  from 7. The %Geometry also defines the boundary conditions for the
 *  electric field calculation. The boundaries that have to be defined
 *  include the six (in 3D) sides of the simulation box. The edges of
 *  the solids are also special kind of boundaries. The boundary
 *  conditions available currently in IBSimu are Dirichlet and Neumann
 *  conditions. Dirichlet boundary condition means that the potential
 *  on the boundary is set to a fixed, known value. Neumann boundary
 *  means that the first derivative of the potential with respect to
 *  the unit outward normal (out of boundary into the simulation
 *  volume) of the surface is set to a preset value. 
 *
 *  The %Geometry class also contains a data mesh, known as the *solid
 *  mesh*, for classifying the nodes into specified types. The solid
 *  mesh, which is internally an array of signed chars, can have
 *  following values and meanings:
 *  
\code
0       = Vacuum
+1...+6 = Dirichlet boundaries at (xmin, xmax, ymin, ymax, zmin, zmax)
-1...-6 = Neumann boundaries at (xmin, xmax, ymin, ymax, zmin, zmax)
+7...   = Solids
-7...   = %Solid edges
\endcode
 *
 *  The solids in the geometry can be defined in several different
 *  ways. The most primitive method of definition is using C++
 *  functions. This kind of solids are known as %FuncSolid. Another
 *  way of definition is through the Constructive %Solid %Geometry
 *  library using the %CSGObjectSolid.
 *
 *  In the following example, a simple 3D geometry is defined with two
 *  solids, the other being a %FuncSolid and the other one a
 *  %CSGObjectSolid, just to demonstrate both. All necessary boundary
 *  conditions and the solid mesh is built: 
\code
#include "geometry.hpp"
#include "func_solid.hpp"
#include "csgobject_solid.hpp"

bool solid7( double x, double y, double z )
{
    return( x*x + y*y + z*z < 0.1*0.1 );
}

CSGObject *solid8( void )
{
    CSGSphere *sphere = new CSGSphere();
    sphere->scale( 0.1, 0.1, 0.1 );
    return( sphere );
}

int main( void )
{
    Geometry geom( MODE_3D, Int3D( 11, 11, 11 ), Vec3D( 0.0, 0.0, 0.0 ), 0.01 );
    Solid *s7 = new FuncSolid( solid7 );
    geom.set_solid( 7, s7 );
    Solid *s8 = new CSGObjectSolid( solid8() );
    geom.set_solid( 7, s8 );
    geom.build_mesh();
    return( 0 );
}
\endcode
 *
 *  See class Geometry for more information in the reference manual.
 */


/*! \brief Boundary condition definition class.
 *
 *  Contains boundary condition type and numerical boundary value.
 *  Dirichlet here means fixed, preset potential at the boundary 
 *  \f[ \phi = \phi_0. \f]
 *  Neumann here means that the first derivative of the potential 
 *  with respect to the unit outward normal (out of solid into the 
 *  vacuum) of the surface is preset 
 *  \f[ - \frac{\partial \phi}{\partial \vec{n}} 
 *      = - \sum_i n_i \frac{\partial \phi}{\partial x_i} 
 *      = q_0. \f]
 */
struct Bound {
    bound_e         type;
    double          val;

    /*! \brief Constructor.
     */
    Bound( bound_e t, double v ) : type(t), val(v) {}

    /*! \brief Constructor for loading boundary condition from a file.
     */
    Bound( std::istream &s ) {
    	type = (bound_e)read_int32( s );
    	val = read_double( s );
    }

    /*! \brief Saves boundary data to stream.
     */
    void save( std::ostream &fout ) const {
	write_int32( fout, type );
	write_double( fout, val );
    }
};


/*! \brief %Geometry defining class.
 *
 *  %Geometry class holds the definitions of the geometry
 *  dimensionality, mesh size and electrode configuration.
 *
 *  Solids are numbered increasingly starting from 7. The solid number
 *  0 is reserved for vacuum and solid numbers from 1 to 6 are
 *  reserved for Dirichlet type boundaries of the bounding box. %Solid
 *  numbers from -1 to -6 are reserved for Neumann type boundaries of
 *  the bounding box. Negative solid numbers starting from -7 are used
 *  for marking electrode edges. The mesh edges are therefore always
 *  'well' defined after running build_mesh().
 *
 *  Bounding box edges are numbered in order xmin, xmax, ymin, ymax,
 *  zmin, xmax.
 */
class Geometry {
    geom_mode_e                _geom_mode; /*!< \brief Geometry mode */
    Int3D                      _size;      /*!< \brief Size of mesh */
    Vec3D                      _origo;     /*!< \brief Location of mesh point (0,0,0) [m] */
    Vec3D                      _max;       /*!< \brief Location of mesh point (size[0]-1,size[1]-1,size[2]-1) [m] */
    double                     _h;         /*!< \brief Length of mesh step [m] */

    int32_t                    _n;         /*!< \brief Number of solids */
    std::vector<const Solid*>  _sdata;     /*!< \brief Array of solid definitions, size \a _n */
    std::vector<Bound>         _bound;     /*!< \brief Array of boundary conditions, size \a _n+6 */

    bool                       _built;     /*!< \brief Is solid mesh array built? */
    signed char               *_smesh;     /*!< \brief Solid mesh array */

public:

    /*! \brief Constructor for geometry class.
     *
     *  Sets geometry mode, mesh cell size \a h, mesh size \a size and
     *  origo \a origo.
     */
    Geometry( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h );

    /*! \brief Constructor for loading geometry from a file.
     */
    Geometry( std::istream &s );

    /*! \brief Destructor for geometry.
     */
    ~Geometry();

    /*! \brief Returns geometry mode.
     */
    geom_mode_e geom_mode( void ) const { return( _geom_mode ); }

    /*! \brief Returns number of dimensions in geometry.
     */
    int32_t dim( void ) const;

    /*! \brief Returns size array of geometry.
     */
    Int3D size( void ) const { return( _size ); }

    /*! \brief Returns size of solid mesh in direction \a i.
     */
    int32_t size( int i ) const { return( _size[i] ); }
   
    /*! \brief Returns number of nodes in the mesh.
     */
    int32_t nodecount( void ) const { return( _size[0]*_size[1]*_size[2] ); }

    /*! \brief Returns origo vector of geometry.
     */
    Vec3D origo( void ) const { return( _origo ); }

    /*! \brief Returns \a i-th component of vector origo.
     */
    double origo( int i ) const { return( _origo[i] ); }

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
    double h( void ) const { return( _h ); }

    /*! \brief Sets solid number \a n to \a s.
     *
     *  Solids have to be defined in incresing order. %Solid number \a
     *  n should be >= 7. This function can also be used to overwrite
     *  a previous solid definition. Pointer to solid \a s is saved to
     *  geometry structure. %Solid will be deleted when geometry is
     *  deleted. The newly defined defined solids default to Dirichlet
     *  boundary with potential zero.
     */
    void set_solid( int32_t n, const Solid *s );

    /*! \brief Returns a const pointer to solid number \a n.
     *
     *  %Solid number \a n should be >= 7.
     */
    const Solid *get_solid( int32_t n ) const;

    /*! \brief Sets boundary condition \a b for solid number \a n.
     *
     *  %Solid number \a n should be > 0 here. Boundary numbers from 1
     *  to 6 are the boundary conditions for the bounding box. Numbers
     *  starting from 7 are the user defined solids. All boundaries of
     *  the simulation box (n <= 6) default to Neumann boundary
     *  condition with derivative value zero. All defined solids (n >=
     *  7) default to Dirichlet boundary with potential zero.
     *
     *  In cylindrical geometry case the \a rmin boundary can be set
     *  to Dirichlet, which means that there is an infinitely thin
     *  wire with a fixed potential at the axis or to Neumann, which
     *  means that the natural boundary for cylindrical axis will be
     *  used.
     */
    void set_boundary( int32_t n, const Bound &b );

    /*! \brief Returns boundary condition for solid number \a n.
     */
    Bound get_boundary( int32_t n ) const;

    /*! \brief Returns a vector of boundary conditions.
     */
    std::vector<Bound> get_boundaries() const;

    /*! \brief Returns 0 if point \a x is vacuum or the number of
     *  solid of \a x is inside a defined solid. Returns a number from
     *  1 to 6 if point \a x is outside the defined geometry.
     */
    int32_t inside( const Vec3D &x ) const;

    /*! \brief Returns true if point \a x is inside solid \a n.
     */
    bool inside( int32_t n, const Vec3D &x ) const;

    /*! \brief Find solid \a n surface location by bracketing.
     *
     *  Searches for the solid \a n surface location on the line
     *  between points \a xin and \a xout by bracketing. Point \a xin
     *  should be inside the solid and point \a xout should be outside
     *  the solid. Function saves the coordinates of the surface to
     *  xsurf and returns parametrical distance (value from 0 to 1)
     *  from xin. %Vector \a xsurf is used as internal work space.
     */
    double bracket_surface( int32_t n, const Vec3D &xin, const Vec3D &xout, Vec3D &xsurf ) const;

    /*! \brief Is the solid mesh built?
     */
    bool built( void ) const { return( _built ); }

    /*! \brief Builds (or rebuilds) the solid mesh from solid
     *  definitions.
     */
    void build_mesh( void );

    /*! \brief Returns a const reference to solid mesh.
     */
    const signed char &mesh( int32_t i ) const { return( _smesh[i] ); }

    /*! \brief Returns a const reference to solid mesh.
     */
    const signed char &mesh( int32_t i, int32_t j ) const {
	return( _smesh[i + j*_size[0]] ); 
    }

    /*! \brief Returns a const reference to solid mesh.
     */
    const signed char &mesh( int32_t i, int32_t j, int32_t k ) const {
	return( _smesh[i + j*_size[0] + k*_size[0]*_size[1]] );
    }

    /*! \brief Returns a reference to solid mesh.
     */
    signed char &mesh( int32_t i ) { return( _smesh[i] ); }

    /*! \brief Returns a reference to solid mesh.
     */
    signed char &mesh( int32_t i, int32_t j ) {
	return( _smesh[i + j*_size[0]] );
    }

    /*! \brief Returns a reference to solid mesh.
     */
    signed char &mesh( int32_t i, int32_t j, int32_t k ) {
	return( _smesh[i + j*_size[0] + k*_size[0]*_size[1]] );
    }

    /*! \brief Returns solid number from solid mesh at \a i, \a j, \a
     *  k or number from 1 to 6 if point is outside mesh.
     */
    signed char mesh_check( int32_t i, int32_t j, int32_t k ) const;

    /*! \brief Saves geometry data to stream.
     */
    void save( std::ostream &s ) const;

    /*! \brief Prints internal data to std::cout.
     */
    void debug_print( void ) const;
};


#endif













