lib_LTLIBRARIES = libibsimu.la

libibsimu_la_LDFLAGS = -version-info 0:1:0
#AM_LDFLAGS = -Wall -O3 -march=athlon64
#AM_CPPFLAGS = -Wall -O3 -march=athlon64

AM_LDFLAGS = -Wall @rdynamic_flag@
AM_CPPFLAGS = -Wall

libibsimu_la_SOURCES = \
bicgstab.cpp \
bicgstab.hpp \
bicgstab_solver.cpp \
bicgstab_solver.hpp \
ccolmatrix.cpp \
ccolmatrix.hpp \
color.hpp \
colormap.cpp \
colormap.hpp \
compmath.cpp \
compmath.hpp \
config.h \
config.hpp \
coordmapper.hpp \
coordmatrix.cpp \
coordmatrix.hpp \
crowmatrix.cpp \
crowmatrix.hpp \
diag_precond.cpp \
diag_precond.hpp \
dxf_solid.cpp \
dxf_solid.hpp \
efield.hpp \
empty_precond.hpp \
epot_efield.cpp \
epot_efield.hpp \
epot_problem.cpp \
epot_problem.hpp \
eqpotgraph.cpp \
eqpotgraph.hpp \
error.cpp \
error.hpp \
fielddiagplot.cpp \
fielddiagplot.hpp \
fielddiagplotter.cpp \
fielddiagplotter.hpp \
fieldgraph.cpp \
fieldgraph.hpp \
fonts.cpp \
fonts.hpp \
frame.cpp \
frame.hpp \
func_solid.cpp \
func_solid.hpp \
geometry.cpp \
geometry.hpp \
geomplot.cpp \
geomplot.hpp \
geomplotter.cpp \
geomplotter.hpp \
graph.hpp \
graph3d.hpp \
gs_solver.cpp \
gs_solver.hpp \
gtkfielddiagdialog.cpp \
gtkfielddiagdialog.hpp \
gtkfielddiagexportdialog.cpp \
gtkfielddiagexportdialog.hpp \
gtkfielddiagwindow.cpp \
gtkfielddiagwindow.hpp \
gtkgeomwindow.cpp \
gtkgeomwindow.hpp \
gtkhardcopy.cpp \
gtkhardcopy.hpp \
gtkparticlediagdialog.cpp \
gtkparticlediagdialog.hpp \
gtkparticlediagwindow.cpp \
gtkparticlediagwindow.hpp \
gtkplotter.cpp \
gtkplotter.hpp \
gtkpreferences.cpp \
gtkpreferences.hpp \
gtkwindow.cpp \
gtkwindow.hpp \
hbio.cpp \
hbio.hpp \
histogram.cpp \
histogram.hpp \
ibsimu.cpp \
ibsimu.hpp \
icon_field_diag.hpp \
icon_hardcopy.hpp \
icon_move.hpp \
icon_particle_diag.hpp \
icon_track.hpp \
icon_zoom_fit.hpp \
icon_zoom_in.hpp \
icon_zoom_out.hpp \
icons.cpp \
icons.hpp \
ilu0_precond.cpp \
ilu0_precond.hpp \
interpolation.cpp \
interpolation.hpp \
label.cpp \
label.hpp \
legend.cpp \
legend.hpp \
lineclip.cpp \
lineclip.hpp \
mat3d.cpp \
mat3d.hpp \
matrix.cpp \
matrix.hpp \
meshgraph.cpp \
meshgraph.hpp \
mvector.cpp \
mvector.hpp \
mydxfblocks.cpp \
mydxfblocks.hpp \
mydxfentities.cpp \
mydxfentities.hpp \
mydxffile.cpp \
mydxffile.hpp \
mydxfheader.cpp \
mydxfheader.hpp \
mydxftables.cpp \
mydxftables.hpp \
palette.cpp \
palette.hpp \
particledatabase.cpp \
particledatabase.hpp \
particlediagplot.cpp \
particlediagplot.hpp \
particlediagplotter.cpp \
particlediagplotter.hpp \
particlegraph.cpp \
particlegraph.hpp \
particleiterator.hpp \
particles.cpp \
particles.hpp \
plotter.cpp \
plotter.hpp \
polysolver.cpp \
polysolver.hpp \
precond.hpp \
problem.hpp \
qrandom.cpp \
qrandom.hpp \
ruler.cpp \
ruler.hpp \
scalarfield.cpp \
scalarfield.hpp \
scharge.cpp \
scharge.hpp \
scheduler.hpp \
solid.hpp \
solidgraph.cpp \
solidgraph.hpp \
solver.hpp \
sort.hpp \
statusprint.cpp \
statusprint.hpp \
symbols.hpp \
timer.cpp \
timer.hpp \
trajectory.cpp \
trajectory.hpp \
trajectorydiagnostics.cpp \
trajectorydiagnostics.hpp \
transformation.cpp \
transformation.hpp \
types.cpp \
types.hpp \
vec3d.cpp \
vec3d.hpp \
vec4d.cpp \
vec4d.hpp \
vectorfield.cpp \
vectorfield.hpp \
xygraph.cpp \
xygraph.hpp

libibsimu_la_includedir=$(includedir)/ibsimu
libibsimu_la_include_HEADERS = \
bicgstab.hpp \
bicgstab_solver.hpp \
ccolmatrix.hpp \
color.hpp \
colormap.hpp \
config.h \
coordmapper.hpp \
coordmatrix.hpp \
crowmatrix.hpp \
diag_precond.hpp \
dxf_solid.hpp \
efield.hpp \
empty_precond.hpp \
epot_efield.hpp \
epot_problem.hpp \
eqpotgraph.hpp \
error.hpp \
fielddiagplot.hpp \
fielddiagplotter.hpp \
fieldgraph.hpp \
file.hpp \
fonts.hpp \
frame.hpp \
func_solid.hpp \
geometry.hpp \
geomplot.hpp \
geomplotter.hpp \
graph.hpp \
graph3d.hpp \
gs_solver.hpp \
gtkfielddiagexportdialog.hpp \
gtkplotter.hpp \
hbio.hpp \
histogram.hpp \
ibsimu.hpp \
ilu0_precond.hpp \
interpolation.hpp \
label.hpp \
legend.hpp \
lineclip.hpp \
mat3d.hpp \
matrix.hpp \
meshgraph.hpp \
mvector.hpp \
mydxfblocks.hpp \
mydxfentities.hpp \
mydxffile.hpp \
mydxfheader.hpp \
mydxftables.hpp \
palette.hpp \
particledatabase.hpp \
particlediagplot.hpp \
particlediagplotter.hpp \
particlegraph.hpp \
particleiterator.hpp \
particles.hpp \
plotter.hpp \
polysolver.hpp \
precond.hpp \
problem.hpp \
qrandom.hpp \
ruler.hpp \
scalarfield.hpp \
scharge.hpp \
scheduler.hpp \
solid.hpp \
solidgraph.hpp \
solver.hpp \
sort.hpp \
statusprint.hpp \
timer.hpp \
trajectory.hpp \
trajectorydiagnostics.hpp \
transformation.hpp \
types.hpp \
vec3d.hpp \
vec4d.hpp \
vectorfield.hpp \
xygraph.hpp

if UMFPACK
libibsimu_la_SOURCES += umfpack_solver.cpp umfpack_solver.hpp
libibsimu_la_include_HEADERS += umfpack_solver.hpp
endif

if CSG
libibsimu_la_SOURCES += csgobject_solid.cpp csgobject_solid.hpp
libibsimu_la_include_HEADERS += csgobject_solid.hpp
endif


CLEANFILES = *~


