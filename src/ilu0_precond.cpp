/*! \file ilu0_precond.cpp
 *  \brief Source code for ilu0_precond.cpp
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

#include "ilu0_precond.hpp"
#include "error.hpp"
#include <stdlib.h>



void ILU0_Precond_constructor( Matrix **_L, Matrix **_U, const CRowMatrix &A )
{
    int i, j, k, pn, qn, rn;
    double mult;
    int l_nz = 0;
    int u_nz = 0;

    // Make checks
    if( A.columns() != A.rows() )
	throw( ErrorDim( ERROR_LOCATION, "matrix not squrare" ) );

    // Calculate number of non-zeros in each matrix and check
    // existance of diagonal.
    int d;
    for( i = 0; i < A.rows(); i++ ) {
	d = 0;
	for( j = A.ptr(i); j < A.ptr(i+1); j++ ) {
	    if( A.col(j) < i )
		l_nz++;
	    else if( A.col(j) > i )
		u_nz++;
	    else {
		u_nz++;
		d++;
	    }
	}
	if( d < 1 )
	    throw( Error( ERROR_LOCATION, "no diagonal entry on row " + to_string(i) ) );
	else if( d > 1 )
	    throw( Error( ERROR_LOCATION, "multiple diagonal entries on row " + to_string(i) ) );
    }

    // Initialize matrices
    CRowMatrix *L = new CRowMatrix( A.rows(), A.rows() );
    CRowMatrix *U = new CRowMatrix( A.rows(), A.rows() );
    L->set_nz( l_nz );
    U->set_nz( u_nz );

    // Build L and U matrices
    L->ptr(0) = U->ptr(0) = 0;
    for( i = 0; i < A.rows(); i++ ) {
	L->ptr(i+1) = L->ptr(i);
	U->ptr(i+1) = U->ptr(i);
    
	for( j = A.ptr(i); j < A.ptr(i+1); j++ )
	    if( A.col(j) < i ) {
		k = L->ptr(i+1)++;
		L->val(k) = A.val(j);
		L->col(k) = A.col(j);
	    } else {
		k = U->ptr(i+1)++;
		U->val(k) = A.val(j);
		U->col(k) = A.col(j);
	    }
    }

    L->order_ascending();
    U->order_ascending();

    // Factor matrix
    for( i = 1; i < A.rows(); i++ ) {
	for( j = L->ptr(i); j < L->ptr(i+1); j++ ) {
	    pn = U->ptr(L->col(j));
	    if( U->val(pn) == 0.0 ) {
		delete L;
		delete U;
		throw( Error( ERROR_LOCATION, "zero pivot in incomplete LU factorization on row " + to_string(i) ) );
	    }
	    mult = (L->val(j) /= U->val(pn));

	    qn = j + 1;
	    rn = U->ptr(i);

	    for( pn++; pn < U->ptr(L->col(j)+1) && U->col(pn) < i; pn++ ) {
		while( qn < L->ptr(i+1) && L->col(qn) < U->col(pn) )
		    qn++;
		if( qn < L->ptr(i+1) && U->col(pn) == L->col(qn) )
		    L->val(qn) -= mult * U->val(pn);
	    }
	    for( ; pn < U->ptr(L->col(j)+1); pn++ ) {
		while( rn < U->ptr(i+1) && U->col(rn) < U->col(pn) )
		    rn++;
		if( rn < U->ptr(i+1) && U->col(pn) == U->col(rn) )
		    U->val(rn) -= mult * U->val(pn);
	    }
	}
    }

    // Save matrices
    *_L = L;
    *_U = U;
}


void ILU0_Precond_constructor( Matrix **_L, Matrix **_U, const CColMatrix &A )
{
    int i, j, k, pn, qn, rn;
    double mult;
    int l_nz = 0;
    int u_nz = 0;

    // Make checks
    if( A.columns() != A.rows() )
	throw( ErrorDim( ERROR_LOCATION, "matrix not squrare" ) );

    // Calculate number of non-zeros in each matrix and check
    // existance of diagonal.
    int d;
    for( i = 0; i < A.rows(); i++ ) {
	d = 0;
	for( j = A.ptr(i); j < A.ptr(i+1); j++ ) {
	    if( A.row(j) > i )
		l_nz++;
	    else if( A.row(j) < i )
		u_nz++;
	    else {
		u_nz++;
		d++;
	    }
	}
	if( d < 1 )
	    throw( Error( ERROR_LOCATION, "no diagonal entry on column " + to_string(i) ) );
	else if( d > 1 )
	    throw( Error( ERROR_LOCATION, "multiple diagonal entries on column " + to_string(i) ) );
    }

    // Initialize matrices
    CColMatrix *L = new CColMatrix( A.rows(), A.rows() );
    CColMatrix *U = new CColMatrix( A.rows(), A.rows() );
    L->set_nz( l_nz );
    U->set_nz( u_nz );

    // Build L and U matrices
    L->ptr(0) = U->ptr(0) = 0;
    for( i = 0; i < A.rows(); i++ ) {
	L->ptr(i+1) = L->ptr(i);
	U->ptr(i+1) = U->ptr(i);
    
	for( j = A.ptr(i); j < A.ptr(i+1); j++ )
	    if( A.row(j) > i ) {
		k = L->ptr(i+1)++;
		L->val(k) = A.val(j);
		L->row(k) = A.row(j);
	    } else {
		k = U->ptr(i+1)++;
		U->val(k) = A.val(j);
		U->row(k) = A.row(j);
	    }
    }

    L->order_ascending();
    U->order_ascending();

    // Factor matrix
    for( i = 0; i < A.rows() - 1; i++ ) {
	mult = U->val(U->ptr(i+1)-1);
	if( mult == 0.0 ) {
	    delete L;
	    delete U;
	    throw( Error( ERROR_LOCATION, "zero pivot in incomplete LU factorization on column " + to_string(i) ) );
	}

	for( j = L->ptr(i); j < L->ptr(i+1); j++ )
	    L->val(j) /= mult;

	for( j = U->ptr(i+1); j < U->ptr(i+2)-1; j++ ) {
	    mult = U->val(j);
	    qn = j + 1;
	    rn = L->ptr(i+1);
	    for( pn = L->ptr(U->row(j)); pn < L->ptr(U->row(j)+1) && L->row(pn) <= i + 1; 
		 pn++ ) {
		while( qn < U->ptr(i+2) && U->row(qn) < L->row(pn) )
		    qn++;
		if( qn < U->ptr(i+2) && L->row(pn) == U->row(qn) )
		    U->val(qn) -= mult * L->val(pn);
	    }
	    for( ; pn < L->ptr(U->row(j)+1); pn++ ) {
		while( rn < L->ptr(i+2) && L->row(rn) < L->row(pn) )
		    rn++;
		if( rn < L->ptr(i+2) && L->row(pn) == L->row(rn) )
		    L->val(rn) -= mult * L->val(pn);
	    }
	}
    }

    // Save matrices
    *_L = L;
    *_U = U;
}


ILU0_Precond::ILU0_Precond( const Matrix &A )
{
    const CRowMatrix *rmat;
    const CColMatrix *cmat;

    if( (rmat = dynamic_cast<const CRowMatrix *>(&A)) != 0 )
	ILU0_Precond_constructor( &_L, &_U, *rmat );
    else if( (cmat = dynamic_cast<const CColMatrix *>(&A)) != 0 )
	ILU0_Precond_constructor( &_L, &_U, *cmat );
    else
	throw( ErrorUnimplemented( ERROR_LOCATION, "ILU0_Precond unimplemented for matrix type" ) );
}


ILU0_Precond::~ILU0_Precond()
{
    delete _L;
    delete _U;
}


void ILU0_Precond::debug_print( void ) const
{
    std::cout << "U = \n" << *_U << "\n\n";
    std::cout << "L = \n" << *_L << "\n\n";
}


void ILU0_Precond::solve( Vector &x, const Vector &b ) const
{
    Vector y(_L->rows());

    if( b.size() != _L->rows() )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    _L->lower_unit_solve( y, b );
    _U->upper_diag_solve( x, y );
}


















