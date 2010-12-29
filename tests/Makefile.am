beam2d_SOURCES = beam2d.cpp
beam3d_SOURCES = beam3d.cpp
beamcyl_SOURCES = beamcyl.cpp
ccolmatrix_SOURCES = ccolmatrix.cpp
colormap_SOURCES = colormap.cpp
crowmatrix_SOURCES = crowmatrix.cpp
efield_SOURCES = efield.cpp
efield_extrpl_SOURCES = efield_extrpl.cpp
efieldtest_2d_SOURCES = efieldtest_2d.cpp
emittance3d_SOURCES = emittance3d.cpp
emittanceconv_SOURCES = emittanceconv.cpp
emittancedef2d_SOURCES = emittancedef2d.cpp
emittancedef3d_SOURCES = emittancedef3d.cpp
emittancedefcyl_SOURCES = emittancedefcyl.cpp
fileplot_SOURCES = fileplot.cpp
fonts_SOURCES = fonts.cpp
fielddiagplot_SOURCES = fielddiagplot.cpp
frame_SOURCES = frame.cpp
frame_extend_range_SOURCES = frame_extend_range.cpp
frame_increase_margin_SOURCES = frame_increase_margin.cpp
geometry_SOURCES = geometry.cpp
gs_solver_SOURCES = gs_solver.cpp
gtkplotter_SOURCES = gtkplotter.cpp
hbio_SOURCES = hbio.cpp
histogram_SOURCES = histogram.cpp
label_SOURCES = label.cpp
lu_test_SOURCES = lu_test.cpp
matvec3dtest_SOURCES = matvec3dtest.cpp
matrixconv_SOURCES = matrixconv.cpp
matrixsolver_SOURCES = matrixsolver.cpp
mvector_SOURCES = mvector.cpp
nsimp_plasmacyl_SOURCES = nsimp_plasmacyl.cpp
palette_SOURCES = palette.cpp
particles_bfield2d_SOURCES = particles_bfield2d.cpp
particles_bfieldcyl_SOURCES = particles_bfieldcyl.cpp
particles_cyl2_SOURCES = particles_cyl2.cpp
particles_cyl_SOURCES = particles_cyl.cpp
particles_diagnostic_SOURCES = particles_diagnostic.cpp
particles_edip_SOURCES = particles_edip.cpp
particles_mirror_SOURCES = particles_mirror.cpp
particles_trajectory_SOURCES = particles_trajectory.cpp
plasma2d_SOURCES = plasma2d.cpp
plasmacyl_SOURCES = plasmacyl.cpp
qrandom_SOURCES = qrandom.cpp
quadrupole_SOURCES = quadrupole.cpp
ruler_SOURCES = ruler.cpp
ruler_inv_SOURCES = ruler_inv.cpp
ruler_tic_bbox_test_SOURCES = ruler_tic_bbox_test.cpp
scharge_out_of_range_SOURCES = scharge_out_of_range.cpp
sclfield_SOURCES = sclfield.cpp
solver1d_edge_SOURCES = solver1d_edge.cpp
solver1d_hard_edge_SOURCES = solver1d_hard_edge.cpp
solver1d_neumann_SOURCES = solver1d_neumann.cpp
solver1d_neumann_order_SOURCES = solver1d_neumann_order.cpp
solver1d_scharge_SOURCES = solver1d_scharge.cpp
solver1d_solid_neumann_conflict_SOURCES = solver1d_solid_neumann_conflict.cpp
solver2d_coax_SOURCES = solver2d_coax.cpp
solver2d_cylindrical_beam_SOURCES = solver2d_cylindrical_beam.cpp
solver2d_smooth_SOURCES = solver2d_smooth.cpp
solver3d_sphere_SOURCES = solver3d_sphere.cpp
solvercyl_beam_SOURCES = solvercyl_beam.cpp
solvercyl_coax_SOURCES = solvercyl_coax.cpp
trajdiaghisto_SOURCES = trajdiaghisto.cpp
vecfield_SOURCES = vecfield.cpp
vecfield_filedef_SOURCES = vecfield_filedef.cpp
vecfield_rotate_SOURCES = vecfield_rotate.cpp
vlasov2d_SOURCES = vlasov2d.cpp

if UMFPACK
umfpack_solver_SOURCES = umfpack_solver.cpp
endif

if CSG
csgobjectsolid_SOURCES = csgobjectsolid.cpp
endif

check_PROGRAMS = \
beam2d \
beam3d \
beamcyl \
ccolmatrix \
colormap \
crowmatrix \
efield \
efield_extrpl \
efieldtest_2d \
emittance3d \
emittanceconv \
emittancedef2d \
emittancedef3d \
emittancedefcyl \
fileplot \
fielddiagplot \
fonts \
frame \
frame_extend_range \
frame_increase_margin \
geometry \
gs_solver \
gtkplotter \
hbio \
histogram \
label \
lu_test \
matvec3dtest \
matrixconv \
matrixsolver \
mvector \
nsimp_plasmacyl \
palette \
particles_bfield2d \
particles_bfieldcyl \
particles_cyl \
particles_cyl2 \
particles_diagnostic \
particles_edip \
particles_mirror \
particles_trajectory \
plasma2d \
plasmacyl \
qrandom \
quadrupole \
ruler \
ruler_inv \
ruler_tic_bbox_test \
scharge_out_of_range \
sclfield \
solver1d_edge \
solver1d_hard_edge \
solver1d_neumann \
solver1d_neumann_order \
solver1d_scharge \
solver1d_solid_neumann_conflict \
solver2d_coax \
solver2d_cylindrical_beam \
solver2d_smooth \
solver3d_sphere \
solvercyl_beam \
solvercyl_coax \
trajdiaghisto \
vecfield \
vecfield_filedef \
vecfield_rotate \
vlasov2d

if UMFPACK
check_PROGRAMS += umfpack_solver
endif

if CSG
check_PROGRAMS += csgobjectsolid 
endif

TESTS = $(check_PROGRAMS)

AM_CPPFLAGS = $(CPPFLAGS) -Wall -g -I$(top_srcdir)/src
AM_LDFLAGS = $(LIBS) -L$(top_srcdir)/src 
# `gsl-config --libs`
LDADD = $(top_builddir)/src/libibsimu.la -lm
if UMFPACK
AM_CPPFLAGS += `pkg-config --cflags umfpack`
AM_LDFLAGS += `pkg-config --libs umfpack`
endif
if CSG
AM_CPPFLAGS += `pkg-config --cflags csg`
AM_LDFLAGS += `pkg-config --libs csg`
endif

EXTRA_DIST = \
beam2d.gpl \
beamcyl.gpl \
efield_extrpl.gpl \
efield.gpl \
matrixsolver.hb \
particles_cyl.gpl \
solver1d_edge.gpl \
solver1d_hard_edge.gpl \
solver1d_neumann.gpl \
solver1d_neumann_order.gpl \
solver1d_scharge.gpl \
solver1d_solid_neumann_conflict.gpl \
solver2d_coax.gpl \
solver2d_cylindrical_beam.gpl \
solvercyl_beam.gpl \
solvercyl_coax.gpl \
vecfield_filedef.in

CLEANFILES = *~ *.png *.eps *.svg *.pdf *.ps *.dat

