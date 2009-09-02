fileplot_SOURCES = fileplot.cpp
particles_edip_SOURCES = particles_edip.cpp
solver1d_neumann_order_SOURCES = solver1d_neumann_order.cpp
fonts_SOURCES = fonts.cpp
label_SOURCES = label.cpp
gtkplotter_SOURCES = gtkplotter.cpp
histogram_SOURCES = histogram.cpp
frame_extend_range_SOURCES = frame_extend_range.cpp
frame_increase_margin_SOURCES = frame_increase_margin.cpp
ruler_SOURCES = ruler.cpp
ruler_tic_bbox_test_SOURCES = ruler_tic_bbox_test.cpp
ruler_inv_SOURCES = ruler_inv.cpp
colormap_SOURCES = colormap.cpp
palette_SOURCES = palette.cpp
beam2d_SOURCES = beam2d.cpp
beam3d_SOURCES = beam3d.cpp
plotgeom_SOURCES = plotgeom.cpp
beamcyl_SOURCES = beamcyl.cpp
ccolmatrix_SOURCES = ccolmatrix.cpp
crowmatrix_SOURCES = crowmatrix.cpp
csgobjectsolid_SOURCES = csgobjectsolid.cpp
efield_SOURCES = efield.cpp
efield_extrpl_SOURCES = efield_extrpl.cpp
efieldtest_2d_SOURCES = efieldtest_2d.cpp
frame_SOURCES = frame.cpp
geometry_SOURCES = geometry.cpp
hbio_SOURCES = hbio.cpp
lu_test_SOURCES = lu_test.cpp
matrixconv_SOURCES = matrixconv.cpp
matrixsolver_SOURCES = matrixsolver.cpp
mvector_SOURCES = mvector.cpp
particles_bfield2d_SOURCES = particles_bfield2d.cpp
particles_bfieldcyl_SOURCES = particles_bfieldcyl.cpp
particles_cyl_SOURCES = particles_cyl.cpp
particles_cyl2_SOURCES = particles_cyl2.cpp
particles_diagnostic_SOURCES = particles_diagnostic.cpp
particles_mirror_SOURCES = particles_mirror.cpp
particles_trajectory_SOURCES = particles_trajectory.cpp
plasma2d_SOURCES = plasma2d.cpp
plasmacyl_SOURCES = plasmacyl.cpp
qrandom_SOURCES = qrandom.cpp
quadrupole_SOURCES = quadrupole.cpp
scharge_out_of_range_SOURCES = scharge_out_of_range.cpp
sclfield_SOURCES = sclfield.cpp
solver1d_edge_SOURCES = solver1d_edge.cpp
solver1d_hard_edge_SOURCES = solver1d_hard_edge.cpp
solver1d_neumann_SOURCES = solver1d_neumann.cpp
solver1d_scharge_SOURCES = solver1d_scharge.cpp
solver2d_coax_SOURCES = solver2d_coax.cpp
solver2d_cylindrical_beam_SOURCES = solver2d_cylindrical_beam.cpp
solver2d_smooth_SOURCES = solver2d_smooth.cpp
solver3d_sphere_SOURCES = solver3d_sphere.cpp
solvercyl_beam_SOURCES = solvercyl_beam.cpp
solvercyl_coax_SOURCES = solvercyl_coax.cpp
solidplot_SOURCES = solidplot.cpp
trajdiaghisto_SOURCES = trajdiaghisto.cpp
vecfield_SOURCES = vecfield.cpp
vecfield_rotate_SOURCES = vecfield_rotate.cpp
vecfield_filedef_SOURCES = vecfield_filedef.cpp
vlasov2d_SOURCES = vlasov2d.cpp

check_PROGRAMS = \
particles_edip \
solver1d_neumann_order \
gtkplotter \
ruler_inv \
ruler_tic_bbox_test \
label \
colormap \
ruler \
palette \
beam2d \
fonts \
beam3d \
frame_extend_range \
frame_increase_margin \
beamcyl \
histogram \
ccolmatrix \
crowmatrix \
csgobjectsolid \
efield \
efield_extrpl \
efieldtest_2d \
frame \
geometry \
hbio \
lu_test \
fileplot \
matrixconv \
matrixsolver \
mvector \
particles_bfield2d \
particles_bfieldcyl \
particles_cyl \
particles_cyl2 \
particles_diagnostic \
particles_mirror \
particles_trajectory \
plasma2d \
plasmacyl \
qrandom \
quadrupole \
scharge_out_of_range \
sclfield \
solver1d_edge \
solver1d_hard_edge \
solver1d_neumann \
solver1d_scharge \
solver2d_coax \
solver2d_cylindrical_beam \
solver2d_smooth \
solver3d_sphere \
solvercyl_beam \
solvercyl_coax \
plotgeom \
solidplot \
trajdiaghisto \
vecfield \
vecfield_rotate \
vecfield_filedef \
vlasov2d

TESTS = $(check_PROGRAMS)

AM_CPPFLAGS = $(CPPFLAGS) -Wall -g -I$(top_srcdir)/src
AM_LDFLAGS = $(LIBS) -I$(top_srcdir)/src `gsl-config --libs`
LDADD = $(top_builddir)/src/libibsimu2.la -lm

EXTRA_DIST = \
beam2d.gpl \
beamcyl.gpl \
efield_extrpl.gpl \
efield.gpl \
particles_cyl.gpl \
solver1d_edge.gpl \
solver1d_hard_edge.gpl \
solver1d_neumann.gpl \
solver1d_scharge.gpl \
solver2d_coax.gpl \
solver2d_cylindrical_beam.gpl \
solvercyl_beam.gpl \
solvercyl_coax.gpl \
matrixsolver.hb \
vecfield_filedef.dat

CLEANFILES = *~ *.png *.ps *.dat 

