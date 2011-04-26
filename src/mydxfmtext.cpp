/*! \file mydxfmtext.cpp
 *  \brief DXF mtext entity
 */

/* Copyright (c) 2010-2011 Taneli Kalvas. All rights reserved.
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


#include "mydxfmtext.hpp"


MyDXFMText::MyDXFMText( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity MTEXT\n";
#endif

    // Default values
    _p[0] = _p[1] = _p[2] = 0.0;
    _text_height = 1.0;
    _rect_width = 1.0;
    _attachment_point = ATTACHMENT_POINT_TOP_LEFT;
    _drawing_direction = DRAWING_DIRECTION_LEFT_TO_RIGHT;

    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 10 )
	    _p[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 )
	    _p[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 )
	    _p[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 40 )
	    _text_height = dxf->group_get_double();
	else if( dxf->group_get_code() == 41 )
	    _rect_width = dxf->group_get_double();

	else if( dxf->group_get_code() == 71 )
	    _attachment_point = dxf->group_get_int16();
	else if( dxf->group_get_code() == 72 )
	    _drawing_direction = dxf->group_get_int16();

	else if( dxf->group_get_code() == 3 || dxf->group_get_code() == 1 )
	    _text += dxf->group_get_string();
	
	else
	    process_group( dxf );
    }

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFMText::explode( MyDXFEntities *ent, class MyDXFFile *dxf, const Transformation *t ) const
{
    MyDXFMText *text = new MyDXFMText( *this );

    // Transform points
    Vec3D p = t->transform_point( text->_p );
    Vec3D h = text->_p + Vec3D(0,text->_text_height,0);
    text->_p = p;
    text->_text_height = norm2( h-p );

    // Add to entities
    ent->add_entity( text );
}


void MyDXFMText::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "MTEXT" );
    write_common( dxf, ostr );

    // Chop into 250 character chunks. Code 1 for last (or only
    // group), others code 3.
    int done = 0;
    int remaining = _text.length();
    while( remaining ) {
	int code = 1;
	if( remaining > 250 )
	    code = 3;
	std::string s =_text.substr( done, 250 );
	done += s.length();
	remaining -= s.length();

	dxf->write_group( code, s.c_str() );
    }

    dxf->write_group( 10, _p[0] );
    dxf->write_group( 20, _p[1] );
    dxf->write_group( 30, _p[2] );

    dxf->write_group( 40, _text_height );
    dxf->write_group( 41, _rect_width );

    dxf->write_group( 71, _attachment_point );
    dxf->write_group( 72, _drawing_direction );

}


void MyDXFMText::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		       const Transformation *t, const double range[4] ) const
{

}


void MyDXFMText::get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const
{
    std::cout << "Warning: bounding box for MText entity not implemented\n";
    min = _p;
    max = _p;

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "MText bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFMText::scale( class MyDXFFile *dxf, double s )
{
    _p *= s;
    _text_height *= s;
    _rect_width *= s;
}


void MyDXFMText::translate( class MyDXFFile *dxf, const Vec3D &dx )
{
    _p += dx;
}


void MyDXFMText::debug_print( std::ostream &os ) const
{
    std::cout << "MTEXT\n";
    std::cout << "  p = {" 
	      << _p[0] << ", " 
	      << _p[1] << ", " 
	      << _p[2] << "}\n";
    std::cout << "  text_height = " << _text_height << "\n"; 
    std::cout << "  rect_width = " << _rect_width << "\n"; 
    std::cout << "  attachment_point = " << _attachment_point << "\n"; 
    std::cout << "  drawing_direction = " << _drawing_direction << "\n"; 
}



