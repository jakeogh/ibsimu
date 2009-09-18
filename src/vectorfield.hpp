/*! \file vectorfield.hpp
 *  \brief Header file for vectorfield.hpp
 */

/* Copyright (c) 2005-2009 Taneli Kalvas. All rights reserved.
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

#ifndef VECTORFIELD_HPP
#define VECTORFIELD_HPP 1


#include "geometry.hpp"
#include "vec3d.hpp"


/*! \brief %Vector field class.
 *
 *  %VectorField is intended to be used as a container for vector
 *  field data. The IBSimu Package uses vector fields for storing
 *  magnetic field. The vector field class provides a subset of vector
 *  operations to copy, sum and scale vector fields. The class also
 *  includes operators for indexed access to mesh elements and
 *  interpolation functions for linearly interpolated smooth field
 *  data. The interpolation function uses linear extrapolation of
 *  field outside the defined mesh. This way it can be ensured that
 *  %VectorField returns sensible values even close by to the edges of
 *  the geometry. 
 *
 *  If the size of mesh is 1 in some direction, then the field is
 *  constant in that direction. Otherwise linear interpolation is
 *  used.
 *
 */
class VectorField {
    geom_mode_e _geom_mode; /*!< \brief Geometry mode */
    Int3D       _size;      /*!< \brief Size of mesh */
    Vec3D       _origo;     /*!< \brief Location of mesh point (0,0,0) [m] */
    double      _h;         /*!< \brief Length of mesh step [m] */
    double      _div_h;     /*!< \brief One over length of mesh step [1/m] */
    double     *_F[3];      /*!< \brief Vector field data in three components
			     *
			     *   If pointer in array is NULL the component 
			     *   is not stored.
			     */

    bool parse_line( const std::string &str, double c[6], double xscale, double fscale, 
		     size_t cdim, size_t fdim, const std::string &filename, size_t linec );

    void transform( int ind[3] );

public:


    /*! \brief Default constructor.
     */
    VectorField() : _geom_mode(MODE_3D), _h(1.0), _div_h(1.0) {
	_F[0] = _F[1] = _F[2] = 0;
    }

    /*! \brief Constructor for geometry from \a g.
     *
     *  Returns a new vector field with geometry parameters (including
     *  mesh size) set from \a g. The field is set to zero in all
     *  locations.
     */
    VectorField( const Geometry &g, bool fout[3] );

    /*! \brief Constructor for set geometry.
     *
     *  Returns a new vector field with geometry set according to
     *  parameters: \a geom_mode is the geometry mode, \a size is the
     *  size of the mesh, \a origo is the location of mesh point
     *  (0,0,0) and \a h is the mesh cell size. The vector field
     *  components marked \a true in array fout are to be defined in
     *  the vector field. Components marked \a false are always
     *  zero. The field is initially set to zero in all locations.
     */
    VectorField( geom_mode_e geom_mode, bool fout[3], Int3D size, 
		 Vec3D origo, double h );

    /*! \brief Constructor for vector field from ascii file.
     *
     *  The vector field for geometry mode \a geom_mode is read in
     *  from file \a filename. The lines starting with # are
     *  skipped. After that the data is read in line-by-line with one
     *  data point per line. The data columns are separated by white
     *  space. The coordinate data are \a (x, y) in 2D, \a (x, r) in
     *  Cyl and \a (x, y, z) in 3D. The field data to be read are
     *  enabled by user with \a fout. The enabled field data
     *  components are read in after the coordinate data from the data
     *  line.
     *
     *  The data points are expected to appear in coordinate sorted
     *  order because the mesh step h is determined from the spatial
     *  difference of first two data points. Spatial coordinates are
     *  multiplied with \a xscale and field components with \a fscale
     *  while read in.
     *
     *  For magnetic fields the particle iterator assumes vector 
     *  field in the following formats:
     *  In 2D:  (x, y, Bz)
     *  In Cyl: (x, r, Bx, Br, Btheta)
     *  In 3D: (x, y, z, Bx, By, Bz)
     *
     */
    VectorField( geom_mode_e geom_mode, bool fout[3], double xscale, 
		 double fscale, std::string filename );

    /*! \brief Copy constructor.
     */
    VectorField( const VectorField &f );

    /*! \brief Destructor.
     */
    ~VectorField();

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

    /*! \brief Returns mesh cell size.
     */
    inline double h( void ) const { return( _h ); }

    /*! \brief Returns inverse mesh cell size.
     */
    inline double div_h( void ) const { return( _div_h ); }

    /*! \brief Translate field in coordinate system.
     */
    void translate( Vec3D x );

    /*! \brief Scale field in coordinate system.
     */
    void scale( double s );

    /*! \brief Rotate field in coordinate system around x-axis.
     */
    void rotate_x( int a );

    /*! \brief Rotate field in coordinate system around y-axis.
     */
    void rotate_y( int a );

    /*! \brief Rotate field in coordinate system around z-axis.
     */
    void rotate_z( int a );

    /*! \brief Clears the field.
     */
    void clear();

    /*! \brief Resets the field geometry.
     *
     *  Sets the field geometry according to the parameters and clears
     *  the field to zero in all locations.
     */
    void reset( geom_mode_e geom_mode, bool fout[3], Int3D size, 
		Vec3D origo, double h );

    /*! \brief Search minimum and maximum vector length values of
     *  vector field.
     */
    void get_minmax( double &min, double &max ) const;

    /*! \brief Get which field components are defined.
     */
    void get_defined_components( bool fout[3] ) const;

    /*! \brief Copy operator.
     */
    VectorField &operator=( const VectorField &f );

    /*! \brief Accumulation operator.
     */
    VectorField &operator+=( const VectorField &f );

    /*! \brief Scaling operator for field values.
     */
    VectorField &operator*=( double x );

    /*! \brief Inverse scaling operator.
     */
    VectorField &operator/=( double x );

    /*! \brief Operator for getting elements of the field mesh.
     */
    const Vec3D operator()( int32_t i ) const;

    /*! \brief Operator for getting elements of the field mesh.
     */
    const Vec3D operator()( int32_t i, int32_t j ) const;

    /*! \brief Operator for getting elements of the field mesh.
     */
    const Vec3D operator()( int32_t i, int32_t j, int32_t k ) const;

    /*! \brief Operator for setting element (i) of the field mesh.
     */
    void set( int32_t i, const Vec3D &v );
    
    /*! \brief Operator for setting element (i,j) of the field mesh.
     */
    void set( int32_t i, int32_t j, const Vec3D &v );

    /*! \brief Operator for setting element (i,j,k) of the field mesh.
     */
    void set( int32_t i, int32_t j, int32_t k, const Vec3D &v );

    /*! \brief Operator for getting linearly interpolated field value
     *  at \a x.
     */
    Vec3D operator()( Vec3D x ) const;

    /*! \brief Prints internal data to std::cout.
     */
    void debug_print( void ) const;
};


#endif













