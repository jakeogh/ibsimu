/*! \file stlfile.hpp
 *  \brief Stereolithography CAD file handling
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

#ifndef STLFILE_HPP
#define STLFILE_HPP 1


#include <fstream>
#include <string>
#include <vector>
#include "vec3d.hpp"


class STLFile {

    class Triangle {

	Vec3D    _normal;
	Vec3D    _p1;
	Vec3D    _p2;
	Vec3D    _p3;
	uint16_t _attr;

	void static read_binary_float_vector( Vec3D &x, std::ifstream &ifstr );
	static void bbox_ppoint( Vec3D &min, Vec3D &max, const Vec3D &p );

    public:

	Triangle( std::ifstream &ifstr );
	~Triangle();

	const Vec3D &normal( void ) const;
	const Vec3D &p1( void ) const;
	const Vec3D &p2( void ) const;
	const Vec3D &p3( void ) const;
	const Vec3D &operator[]( int i ) const;
	
	void update_bbox( Vec3D &min, Vec3D &max ) const;

	int ray_cross( const Vec3D &x, const Vec3D &l ) const;

	void debug_print( std::ostream &os ) const;
    };

    class VTriangle {

	uint32_t _v[3];
	Vec3D    _normal;

    public:

	VTriangle( uint32_t v1, uint32_t v2, uint32_t v3, const Vec3D &normal );
	VTriangle( const uint32_t v[3], const Vec3D &normal );
	~VTriangle();

	const Vec3D &normal( void ) const;
	const uint32_t &operator[]( int i ) const;

	bool inside( const Vec3D &x ) const;

	void debug_print( std::ostream &os ) const;

	/*! \brief Outputting to stream.
	 */
	friend std::ostream &operator<<( std::ostream &os, const VTriangle &vtri );
    };

    bool                   _ascii;
    std::vector<Triangle>  _triangle;  // Original triangle data

    std::vector<Vec3D>     _vertex;    // Vertex list
    std::vector<VTriangle> _vtri;      // Vertex made triangles

    void read_binary( std::ifstream &ifstr );
    void build_vtriangle_data( void );

public:

    STLFile( const std::string &filename );

    ~STLFile();

    size_t size( void );

    bool inside( const Vec3D &x, double eps = 1.0e-6 ) const;

    /*! \brief Return bounding box of entity
     */
    void get_bbox( Vec3D &min, Vec3D &max ) const;

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;

    /*! \brief Outputting to stream.
     */
    friend std::ostream &operator<<( std::ostream &os, const VTriangle &vtri );
};


#endif


