/*! \file geometry.hpp
 *  \brief %Geometry definition
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP 1


#include <stdint.h>
#include <vector>
#include <iostream>
#include "file.hpp"
#include "vec3d.hpp"
#include "solid.hpp"
#include "mesh.hpp"
#include "types.hpp"


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
struct Bound 
{
    bound_e         type;
    double          val;

    /*! \brief Constructor.
     */
    Bound( bound_e t, double v ) : type(t), val(v) {}

    /*! \brief Constructor for loading boundary condition from a file.
     */
    Bound( std::istream &is ) {
    	type = (bound_e)read_int32( is );
    	val = read_double( is );
    }

    /*! \brief Saves data to stream \a os.
     */
    void save( std::ostream &os ) const {
	write_int32( os, type );
	write_double( os, val );
    }

    /*! \brief Outputting to stream.
     */
    friend std::ostream &operator<<( std::ostream &os, const Bound &b );
};


/*
#define SMESH_NODE_ID_MASK             0xC0000000 // 1100...
#define SMESH_NODE_ID_NEAR_SOLID       0x00000000 // 0000...
#define SMESH_NODE_ID_PURE_VACUUM      0x40000000 // 0100...
#define SMESH_NODE_ID_NEUMANN          0x80000000 // 1000...
#define SMESH_NODE_ID_DIRICHLET        0xC0000000 // 1100...
*/


#define SMESH_NODE_ID_MASK             0xE0000000 // 111...

#define SMESH_NODE_ID_PURE_VACUUM      0x00000000 // 000...
#define SMESH_NODE_ID_NEAR_SOLID       0x20000000 // 001...
#define SMESH_NODE_ID_NEUMANN          0x40000000 // 010...
#define SMESH_NODE_ID_ROUGH_BOUNDARY   0x60000000 // 011...

#define SMESH_NODE_ID_PURE_VACUUM_FIX  0x80000000 // 100...
#define SMESH_NODE_ID_NEAR_SOLID_FIX   0xA0000000 // 101...
#define SMESH_NODE_ID_DIRICHLET        0xC0000000 // 110...
#define SMESH_NODE_ID_FINE_BOUNDARY    0xE0000000 // 111...

#define SMESH_NODE_FIXED               0x80000000 // 100...


#define SMESH_BOUNDARY_NUMBER_MASK     0x000000FF // limit to 0-255
#define SMESH_NEAR_SOLID_INDEX_MASK    0x1FFFFFFF // limit to 0-2^29 (5.4e8)


/*! \brief %Geometry defining class.
 *
 *  %Geometry class holds the definitions of the geometry
 *  dimensionality, mesh size and electrode configuration. Also it
 *  contains a signed char array for information about the type of
 *  each node. This array is known as the solid mesh.
 *
 *  The integer numbers in the solid mesh have the following meanings:
 *  The solid number 0 is reserved for vacuum and solid numbers from 1
 *  to 6 are reserved for Dirichlet type boundaries of the bounding
 *  box. %Solid numbers from -1 to -6 are reserved for Neumann type
 *  boundaries of the bounding box. Negative solid numbers starting
 *  from -7 are used for marking electrode edges and the positive
 *  numbers starting from 7 are used to mark the interior points of
 *  the electrodes.
 *
 *  The mesh nodes are marked using the following logic: First nodes,
 *  which are inside electrodes are marked solid (>=7). If a point is
 *  inside several solids, the highest solid number is marked. Other
 *  points are left as vacuum nodes (0). As the next step the solid
 *  nodes are mapped to find nodes which have vacuum as closest
 *  neighbour along any of the axes (not diagonal). These nodes are
 *  marked as solid edges (<=-7). As the last step, the vacuum nodes
 *  at the simulation box boundary are marked either as Neumann (<0
 *  and >-7) or Dirichlet (>0 and <7).
 *
 *  Starting from 1.0.3: A. Mark solids, B. Mark Neumann and Dirichlet
 *  boundaries, C. Mark edges taking in account that Neumann = Vacuum
 *  and Dirichlet != Vacuum.
 *
 *  Bounding box edges are numbered in order xmin, xmax, ymin, ymax,
 *  zmin, xmax.
 */
class Geometry : public Mesh
{
    uint32_t                   _n;         /*!< \brief Number of solids */
    std::vector<const Solid*>  _sdata;     /*!< \brief Array of solid definitions, size \a _n */
    std::vector<Bound>         _bound;     /*!< \brief Array of boundary conditions, size \a _n+6 */

    bool                       _built;     /*!< \brief Is solid mesh array built? */
    uint32_t                  *_smesh;     /*!< \brief Solid mesh array. */
    std::vector<uint8_t>       _nearsolid; /*!< \brief Near solid data. */

    
    /*! \brief Check if node is solid (n>=7).
     *
     *  Returns 0 if node is not solid (vacuum, outside mesh or
     *  simulation box boundary). If node is solid, the solid number
     *  >= 7 is returned.
     */
    uint32_t is_solid( int32_t i, int32_t j, int32_t k ) const;

    /*! \brief Add a near solid distance to near solid data.
     *
     *  Subroutine of add_near_solid_entry().
     */
    void add_near_solid_distance( std::vector<uint8_t> &ndist, uint8_t dist );

    /*! \brief Add an entry to near solid data.
     *
     *  Appends a new entry to near solid data and updated \a
     *  near_solid_index to point to the next byte after the end of
     *  data.
     */
    void add_near_solid_entry( uint32_t &near_solid_index, int32_t i, int32_t j, int32_t k );

    /*! \brief Bracket solid surface between node points
     *
     *  Finds surface of solid \a solid between the outer node (\a i,
     *  \a j, \a k) and the inner node in \a sign direction (+/-1) of
     *  of coordinate axis \a coord (0,1 or 2). Returns parametric
     *  distance (between 0 and 255) from the outer node to the inner
     *  node.
     */
    uint8_t bracket_ndist( int32_t i, int32_t j, int32_t k, int32_t solid, int sign, int coord ) const;

    /*! \brief Check mesh definition validity.
     */
    void check_definition();

    Vec3D surface_normal_2d( const Vec3D &x ) const;
    Vec3D surface_normal_3d( const Vec3D &x ) const;

public:

    /*! \brief Constructor for geometry class.
     *
     *  Sets geometry mode, mesh cell size \a h, mesh size \a size and
     *  origo \a origo.
     */
    Geometry( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h );

    /*! \brief Constructor for loading geometry from a stream \a is.
     */
    Geometry( std::istream &is );

    /*! \brief Destructor for geometry.
     */
    ~Geometry();

    /*! \brief Return number of solids.
     */
    uint32_t number_of_solids() const;

    /*! \brief Return number of boundaries.
     *
     *  Always returns at least \a n >= 6.
     */
    uint32_t number_of_boundaries() const;
    
    /*! \brief Sets solid number \a n to \a s.
     *
     *  Solids have to be defined in incresing order. %Solid number \a
     *  n should be >= 7. This function can also be used to overwrite
     *  a previous solid definition. Pointer to solid \a s is saved to
     *  geometry structure. %Solid will not be deleted when geometry is
     *  deleted. The newly defined defined solids default to Dirichlet
     *  boundary with potential zero.
     */
    void set_solid( uint32_t n, const Solid *s );

    /*! \brief Returns a const pointer to solid number \a n.
     *
     *  %Solid number \a n should be >= 7.
     */
    const Solid *get_solid( uint32_t n ) const;

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
    void set_boundary( uint32_t n, const Bound &b );

    /*! \brief Returns boundary condition for solid number \a n.
     */
    Bound get_boundary( uint32_t n ) const;

    /*! \brief Returns a vector of boundary conditions.
     */
    std::vector<Bound> get_boundaries() const;

    /*! \brief Returns 0 if point \a x is vacuum or the number of
     *  solid of \a x is inside a defined solid. Returns a number from
     *  1 to 6 if point \a x is outside the defined geometry. If the
     *  point is inside several defined solids, the solid with the
     *  highest solid number is returned.
     */
    uint32_t inside( const Vec3D &x ) const;

    /*! \brief Returns true if point \a x is inside solid \a n.
     */
    bool inside( uint32_t n, const Vec3D &x ) const;

    /*! \brief Find solid \a n surface location by bracketing.
     *
     *  Searches for the solid \a n surface location on the line
     *  between points \a xin and \a xout by bracketing. Point \a xin
     *  should be inside the solid and point \a xout should be outside
     *  the solid. Function saves the coordinates of the surface to
     *  xsurf and returns parametrical distance (value from 0 to 1)
     *  from xin.
     */
    double bracket_surface( uint32_t n, const Vec3D &xin, const Vec3D &xout, Vec3D &xsurf ) const;

    /*! \brief Find surface outward normal at location \a x.
     *
     *  Returns zero vector on failure.
     */
    Vec3D surface_normal( const Vec3D &x ) const;

    /*! \brief Is the solid mesh built?
     */
    bool built( void ) const { return( _built ); }

    /*! \brief Builds (or rebuilds) the solid mesh from solid
     *  definitions.
     */
    void build_mesh( void );

    /*! \brief Returns a const reference to solid mesh array.
     */
    const uint32_t &mesh( int32_t i ) const { return( _smesh[i] ); }

    /*! \brief Returns a const reference to solid mesh array.
     */
    const uint32_t &mesh( int32_t i, int32_t j ) const {
	return( _smesh[i + j*_size[0]] ); 
    }

    /*! \brief Returns a const reference to solid mesh array.
     */
    const uint32_t &mesh( int32_t i, int32_t j, int32_t k ) const {
	return( _smesh[i + j*_size[0] + k*_size[0]*_size[1]] );
    }

    /*! \brief Returns a reference to solid mesh array.
     */
    uint32_t &mesh( int32_t i ) { return( _smesh[i] ); }

    /*! \brief Returns a reference to solid mesh array.
     */
    uint32_t &mesh( int32_t i, int32_t j ) {
	return( _smesh[i + j*_size[0]] );
    }

    /*! \brief Returns a reference to solid mesh array.
     */
    uint32_t &mesh( int32_t i, int32_t j, int32_t k ) {
	return( _smesh[i + j*_size[0] + k*_size[0]*_size[1]] );
    }

    /*! \brief Returns number from solid mesh array.
     *
     *  Returns number from solid mesh array at \a i, \a j, \a
     *  k or Dirichlet boundary number (1-6) if point is outside mesh.
     */
    uint32_t mesh_check( int32_t i, int32_t j, int32_t k ) const;

    /*! \brief Returns true if node is a potential near solid point.
     *
     *  Returns true if any of the neighbouring points is a solid
     *  point (Dirichlet with solid number >= 7).
     */
    bool is_near_solid( int32_t i, int32_t j, int32_t k ) const;

    /*! \brief Returns a const pointer to start of near solid data for
     *  node (\a i, \a j, \a k). The first byte contains the bit flags
     *  for the existance of neighbouring solids. From bit 0 to bit 5
     *  the boolean flags are for directions: xmin, xmax, ymin, ymax,
     *  zmin, zmax. The next bytes contain the parametric distances of
     *  the solid surfaces from the node in each direction. Only the
     *  directions with set bit flag are saved to data. The distances
     *  are saved in the same order as the flags (from xmin to
     *  zmax). The distance information is an unsigned 8-bit integer
     *  (0 to 255), where 0 means distance 0.0 and 255 means 1.0.
     */
    /*
    const uint8_t *nearsolid_ptr( int32_t i, int32_t j, int32_t k ) const {
	return( &_nearsolid[_smesh[i + j*_size[0] + k*_size[0]*_size[1]] & 
			    SMESH_NEAR_SOLID_INDEX_MASK] );
    }
    */
    const uint8_t *nearsolid_ptr( int32_t index ) const {
	return( &_nearsolid[index] );
    }

    /*! \brief Returns distance of solid boundary from point.
     *
     *  Returns the distance (0 to 255) of near solid into direction
     *  \a dir from near solid point at (\a i, \a j, \a k). The
     *  direction \a dir is an integer from 0 to 5, with 0 meaning -x,
     *  1 meaning +x, 2 meaning -y, 3 meaning +y, 4 meaning -z and 5
     *  meaning +z. If the node at (\a i, \a j, \a k) is not a near
     *  solid node or if there is no solid nearby in the direction an
     *  error will be thrown.
     */
    uint8_t solid_dist( uint32_t i, uint32_t j, uint32_t k, uint32_t dir ) const;

    /*! \brief Returns distance of solid boundary from point.
     *
     *  Same as solid_dist() above, just using one dimensional index
     *  for mesh.
     */
    uint8_t solid_dist( uint32_t i, uint32_t dir ) const;

    /*! \brief Saves data to a new file \a filename.
     */
    void save( const std::string &filename ) const;

    /*! \brief Saves data to stream \a os.
     */
    void save( std::ostream &os ) const;

    /*! \brief Print debugging information to stream \a os.
     */
    void debug_print( std::ostream &os ) const;
};


#endif

