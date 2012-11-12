/*! \file stlfile.hpp
 *  \brief Stereolithography CAD file handling
 */

/* Copyright (c) 2011-2012 Taneli Kalvas. All rights reserved.
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

#ifndef STLFILE_HPP
#define STLFILE_HPP 1


#include <fstream>
#include <string>
#include <vector>
#include <stdint.h>
#include "vtriangle.hpp"
#include "vec3d.hpp"


class STLFile {

    class Triangle {

	Vec3D    _normal;
	Vec3D    _p[3];
	uint16_t _attr;

	void static read_binary_float_vector( Vec3D &x, std::ifstream &ifstr );
	void static read_ascii_float_vector( Vec3D &x, const char *buf, 
					     const std::string &filename, int linec );
	static void bbox_ppoint( Vec3D &min, Vec3D &max, const Vec3D &p );

    public:

	Triangle( std::ifstream &ifstr );
	Triangle( std::ifstream &ifstr, const char *buf, const std::string &filename, int &linec );
	~Triangle();

	uint16_t attr( void ) const;
	const Vec3D &normal( void ) const;
	const Vec3D &operator[]( int i ) const {
	    return( _p[i] );
	}
	
	void update_bbox( Vec3D &min, Vec3D &max ) const;
	void debug_print( std::ostream &os ) const;
    };

    std::string            _filename;
    bool                   _ascii;
    std::vector<Triangle>  _tri;       // Original triangles
    std::vector<Vec3D>     _vertex;    // Vertices
    std::vector<VTriangle> _vtri;      // Vertex made triangles

    Vec3D                  _offset;    // Amount of offset on vertices
    Vec3D                  _bbox[2];   // bbox: min, max

    double                 _vertex_matching_eps;
    double                 _signed_volume_eps;

    std::vector<char>      _vpos;
    std::vector<char>      _vneg;

    void read_binary( std::ifstream &ifstr );
    void read_ascii( std::ifstream &ifstr );
    void build_vtriangle_data( void );
    void check_vtriangle_data( void );
    void offset_vtriangle_data( void );
    int classify_original_tetrahedron( int ss, const Vec3D &p, 
				       const Vec3D &q1, const Vec3D &q2, const Vec3D &q3 );
    int signvol4( const Vec3D &q0, const Vec3D &q1, 
		  const Vec3D &q2, const Vec3D &q3 );
    int signvol3( const Vec3D &q1, const Vec3D &q2, const Vec3D &q3 );

    void calculate_bbox( void );

public:

    /*! \brief Constructor for STL data from file.
     *
     *  Reads either binary or ascii STL-file from \a filename.
     *  Triangle vertices are connected if closer that \a
     *  vertex_matching_eps together (absolute distance in
     *  meters). During inside() evaluation a tetrahedron volume less
     *  than \a signed_volume_eps is judged to be in the limits of
     *  numerical accuracy.
     */
    STLFile( const std::string &filename, 
	     double vertex_matching_eps = 1.0e-9, 
	     double signed_volume_eps = 1.0e-15 );

    /*! \brief Constructor for $STLFile from triangle and vertex data
     */
    STLFile( const std::vector<Vec3D> &vertex,
	     const std::vector<VTriangle> &triangle );

    /*! \brief Destructor.
     */
    ~STLFile();

    /*! \brief Write to file.
     */
    void save( const std::string &filename, bool ascii = true ) const;

    /*! \brief Return number of vertices.
     */
    size_t vertexc( void );

    /*! \brief Return vertex \a i coordinates.
     */
    Vec3D vertex( uint32_t i ) const;

    /*! \brief Return number of triangles.
     */
    size_t trianglec( void );

    /*! \brief Return vertex triangle \a i.
     */
    const VTriangle &vtriangle( uint32_t i ) const;

    /*! \brief Return if point \a x is inside solid.
     *
     *  Uses algorithm from R. Segura, F. R. Feito, J. Ruiz de Miras,
     *  C. Ogayar and J. C. Torres, "An Efficient Point Classification
     *  Algorithm for Triangle Meshes", Journal of Graphics, GPU, and
     *  Game Tools 10, issue 3, 2005.
     */
    bool inside( const Vec3D &x );

    /*! \brief Return bounding box of entity
     */
    void get_bbox( Vec3D &min, Vec3D &max ) const;

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;
};


#endif


