#include "fielddiagplot.hpp"


FieldDiagPlot::FieldDiagPlot( Frame *frame, const Geometry *geom )
    : _frame(frame), _geom(geom), _epot(NULL), _efield(NULL), _scharge(NULL), _bfield(NULL), 
      _N(100)
{
    _diag[0] = FIELDD_DIAG_EPOT;
    _diag[1] = FIELDD_DIAG_NONE;
    _loc[0] = FIELDD_LOC_DIST;
    _loc[1] = FIELDD_LOC_NONE;

    _graph[0] = NULL;
    _graph[1] = NULL;
}


FieldDiagPlot::~FieldDiagPlot()
{
    if( _graph[0] )
	delete _graph[0];
    if( _graph[1] )
	delete _graph[1];
}


void FieldDiagPlot::build_plot( void )
{
    // Build coordinate points
    std::vector<double> coord[4];
    coord[0].reserve( _N );
    coord[1].reserve( _N );
    coord[2].reserve( _N );
    coord[3].reserve( _N );
    for( size_t a = 0; a < 3; a++ ) {
	for( size_t b = 0; b < _N; b++ ) {
	    coord[a].push_back( _x1[a] + (b/(_N-1.0))*(_x2[a]-_x1[a]) );
	}
    }
    for( size_t b = 0; b < _N; b++ ) {
	Vec3D d = (b/(_N-1.0))*(_x2-_x1);
	coord[3].push_back( d.norm2() );
    }

    // Set x-ranges manually
    _frame->ruler_autorange_enable( PLOT_AXIS_X1, false, false );
    _frame->ruler_autorange_enable( PLOT_AXIS_X2, false, false );

    //_frame->ruler_autorange_enable( PLOT_AXIS_Y1, false, false );
    //_frame->ruler_autorange_enable( PLOT_AXIS_Y2, false, false );

    // Enable y-axes according to use
    if( _diag[0] != FIELDD_DIAG_NONE )
	_frame->force_enable_ruler( PLOT_AXIS_Y1, true );
    if( _diag[1] != FIELDD_DIAG_NONE )
	_frame->force_enable_ruler( PLOT_AXIS_Y2, true );

    // X data
    PlotAxis xaxis_use;
    std::vector<double> *datax = NULL;
    for( size_t b = 0; b < 2; b++ ) {

	// Check if xaxis disabled
	if( _loc[b] == FIELDD_LOC_NONE ) {
	    //std::cout << "Disabled\n";
	    continue;
	}
	
	// Select xaxis accordingly
	PlotAxis xaxis;
	if( b == 0 )
	    xaxis = PLOT_AXIS_X1;
	else
	    xaxis = PLOT_AXIS_X2;

	// Force enable ruler
	_frame->force_enable_ruler( xaxis, true );

	// Get datax pointer and set xaxis properties
	std::vector<double> *dataxt = NULL;
	if( _loc[b] == FIELDD_LOC_DIST ) { // Distance
	    _frame->set_axis_label( xaxis, "Distance (m)" );
	    _frame->set_ranges( xaxis, 0.0, coord[3][_N-1] );
	    dataxt = &coord[3];
	    //std::cout << "Distance: " << 0.0 << " to " << coord[3][_N-1] << "\n";
	} else if( _loc[b] == FIELDD_LOC_X ) { // X
	    _frame->set_axis_label( xaxis, "x (m)" );
	    _frame->set_ranges( xaxis, _x1[0], _x2[0] );
	    dataxt = &coord[0];
	    //std::cout << "X: " << _x1[0] << " to " << _x2[0] << "\n";
	} else if( _loc[b] == FIELDD_LOC_Y ) { // Y
	    if( _geom->geom_mode() == MODE_CYL )
		_frame->set_axis_label( xaxis, "r (m)" );
	    else
		_frame->set_axis_label( xaxis, "y (m)" );
	    _frame->set_ranges( xaxis, _x1[1], _x2[1] );
	    dataxt = &coord[1];
	    //std::cout << "Y: " << _x1[1] << " to " << _x2[1] << "\n";
	} else if( _loc[b] == FIELDD_LOC_Z ) { // Z
	    _frame->set_axis_label( xaxis, "z (m)" );
	    _frame->set_ranges( xaxis, _x1[2], _x2[2] );
	    dataxt = &coord[2];
	    //std::cout << "Z: " << _x1[2] << " to " << _x2[2] << "\n";
	}

	// Use dataxt if values not constant
	if( dataxt != NULL && (*dataxt)[0] != (*dataxt)[_N-1] ) {
	    datax = dataxt;
	    xaxis_use = xaxis;
	}
    }
    if( datax == NULL )
	throw( Error( ERROR_LOCATION, "both x-axes undefined" ) );


    // Set data
    for( size_t a = 0; a < 2; a++ ) {

	// Check if there is data for yaxis
	if( _diag[a] == 0 )
	    continue;

	// Select yaxis accordingly
	PlotAxis yaxis;
	if( a == 0 )
	    yaxis = PLOT_AXIS_Y1;
	else
	    yaxis = PLOT_AXIS_Y2;

	// Set yaxis properties and data
	std::vector<double> datay;
	datay.reserve( _N );
	if( _diag[a] == FIELDD_DIAG_EPOT ) {
	    _frame->set_axis_label( yaxis, "\\phi  (V)" );
	    if( _epot == NULL )
		throw( Error( ERROR_LOCATION, "epot undefined" ) );		
	    for( size_t b = 0; b < _N; b++ ) {
		datay.push_back( (*_epot)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) ) );
	    }

	} else if( _diag[a] == FIELDD_DIAG_EFIELD ) {
	    _frame->set_axis_label( yaxis, "|E| (V/m)" );
	    if( _efield == NULL )
		throw( Error( ERROR_LOCATION, "efield undefined" ) );		
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D E = (*_efield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		datay.push_back( E.norm2() );
	    }

	} else if( _diag[a] == FIELDD_DIAG_EFIELD_X ) {
	    _frame->set_axis_label( yaxis, "E_x (V/m)" );
	    if( _efield == NULL )
		throw( Error( ERROR_LOCATION, "efield undefined" ) );		
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D E = (*_efield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		datay.push_back( E[0] );
	    }

	} else if( _diag[a] == FIELDD_DIAG_EFIELD_Y ) {
	    if( _geom->geom_mode() == MODE_CYL )
		_frame->set_axis_label( yaxis, "E_r (V/m)" );
	    else
		_frame->set_axis_label( yaxis, "E_y (V/m)" );
	    if( _efield == NULL )
		throw( Error( ERROR_LOCATION, "efield undefined" ) );		
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D E = (*_efield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		datay.push_back( E[1] );
	    }

	} else if( _diag[a] == FIELDD_DIAG_EFIELD_Z ) {
	    _frame->set_axis_label( yaxis, "E_z (V/m)" );
	    if( _efield == NULL )
		throw( Error( ERROR_LOCATION, "efield undefined" ) );		
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D E = (*_efield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		datay.push_back( E[2] );
	    }

	} else if( _diag[a] == FIELDD_DIAG_SCHARGE ) {
	    _frame->set_axis_label( yaxis, "\\rho  (C/m^3)" );
	    if( _scharge == NULL )
		throw( Error( ERROR_LOCATION, "scharge undefined" ) );
	    for( size_t b = 0; b < _N; b++ ) {
		datay.push_back( (*_scharge)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) ) );
	    }

	} else if( _diag[a] == FIELDD_DIAG_BFIELD ) {
	    _frame->set_axis_label( yaxis, "B (T)" );
	    if( _bfield == NULL )
		throw( Error( ERROR_LOCATION, "bfield undefined" ) );
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D B = (*_bfield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		datay.push_back( B.norm2() );
	    }

	} else if( _diag[a] == FIELDD_DIAG_BFIELD_X ) {
	    _frame->set_axis_label( yaxis, "B_x (T)" );
	    if( _bfield == NULL )
		throw( Error( ERROR_LOCATION, "bfield undefined" ) );
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D B = (*_bfield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		datay.push_back( B[0] );
	    }

	} else if( _diag[a] == FIELDD_DIAG_BFIELD_Y ) {
	    if( _geom->geom_mode() == MODE_CYL )
		_frame->set_axis_label( yaxis, "B_r (T)" );
	    else
		_frame->set_axis_label( yaxis, "B_y (T)" );
	    if( _bfield == NULL )
		throw( Error( ERROR_LOCATION, "bfield undefined" ) );
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D B = (*_bfield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		datay.push_back( B[1] );
	    }

	} else if( _diag[a] == FIELDD_DIAG_BFIELD_Z ) {
	    if( _geom->geom_mode() == MODE_CYL )
		_frame->set_axis_label( yaxis, "B_\\theta  (T)" );
	    else
		_frame->set_axis_label( yaxis, "B_z (T)" );
	    if( _bfield == NULL )
		throw( Error( ERROR_LOCATION, "bfield undefined" ) );
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D B = (*_bfield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		datay.push_back( B[2] );
	    }
	}

	// Make graph
	_graph[a] = new XYGraph( *datax, datay );
	_graph[a]->set_line_style( XYGRAPH_LINE_SOLID );
	_graph[a]->set_point_style( XYGRAPH_POINT_DISABLE, true, 0.0 );
	if( a == 0 )
	    _graph[a]->set_color( Color(1,0,0) );
	else
	    _graph[a]->set_color( Color(0,0,1) );
	_frame->add_graph( xaxis_use, yaxis, _graph[a] );
    }

    /*
    // Save default "zoom fit" ranges
    _frame->get_ranges( PLOT_AXIS_X1, _range[0], _range[2] );
    _frame->get_ranges( PLOT_AXIS_Y1, _range[1], _range[3] );
    _frame->get_ranges( PLOT_AXIS_X2, _range[4], _range[6] );
    _frame->get_ranges( PLOT_AXIS_Y2, _range[5], _range[7] );
    */

}
