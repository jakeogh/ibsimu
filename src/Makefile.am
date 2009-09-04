lib_LTLIBRARIES = libibsimu.la

libibsimu_la_LDFLAGS = -version-info 0:1:0
#AM_LDFLAGS = -Wall -O3 -march=athlon64
#AM_CPPFLAGS = -Wall -O3 -march=athlon64

AM_LDFLAGS = -Wall -Werror
AM_CPPFLAGS = -Wall -Werror

libibsimu_la_SOURCES = \
fielddiagplot.hpp \
fielddiagplot.cpp \
fielddiagplotter.hpp \
fielddiagplotter.cpp \
particlediagplot.hpp \
particlediagplot.cpp \
particlediagplotter.hpp \
particlediagplotter.cpp \
geomplot.hpp \
geomplot.cpp \
geomplotter.hpp \
geomplotter.cpp \
plotter.hpp \
plotter.cpp \
statusprint.hpp \
statusprint.cpp \
fieldgraph.hpp \
fieldgraph.cpp \
fonts.hpp \
fonts.cpp \
interpolation.hpp \
interpolation.cpp \
gtkpreferences.cpp \
gtkpreferences.hpp \
gtkgeomwindow.cpp \
gtkgeomwindow.hpp \
gtkwindow.cpp \
gtkwindow.hpp \
gtkhardcopy.cpp \
gtkhardcopy.hpp \
gtkplotter.cpp \
gtkplotter.hpp \
meshgraph.cpp \
meshgraph.hpp \
ruler.cpp \
ruler.hpp \
histogram.hpp \
histogram.cpp \
legend.cpp \
legend.hpp \
palette.cpp \
palette.hpp \
colormap.cpp \
colormap.hpp \
bicgstab.cpp \
bicgstab.hpp \
bicgstab_solver.cpp \
bicgstab_solver.hpp \
color.hpp \
ccolmatrix.cpp \
ccolmatrix.hpp \
config.hpp \
label.cpp \
label.hpp \
coordmatrix.cpp \
coordmatrix.hpp \
crowmatrix.cpp \
crowmatrix.hpp \
csgobject_solid.cpp \
csgobject_solid.hpp \
diag_precond.cpp \
diag_precond.hpp \
efield.hpp \
types.cpp \
types.hpp \
empty_precond.hpp \
epot_efield.cpp \
epot_efield.hpp \
epot_problem.cpp \
epot_problem.hpp \
eqpotgraph.cpp \
eqpotgraph.hpp \
icons.cpp \
icons.hpp \
error.hpp \
frame.cpp \
frame.hpp \
graph.hpp \
graph3d.hpp \
coordmapper.hpp \
func_solid.cpp \
func_solid.hpp \
gauss_seidel.cpp \
gauss_seidel.hpp \
geometry.cpp \
geometry.hpp \
hbio.cpp \
hbio.hpp \
ilu0_precond.cpp \
ilu0_precond.hpp \
lineclip.cpp \
lineclip.hpp \
matrix.cpp \
matrix.hpp \
mvector.cpp \
mvector.hpp \
particledatabase.cpp \
particledatabase.hpp \
particleiterator.hpp \
particles.cpp \
particles.hpp \
particlegraph.cpp \
particlegraph.hpp \
polysolver.cpp \
polysolver.hpp \
precond.hpp \
problem.hpp \
qrandom.cpp \
qrandom.hpp \
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
timer.cpp \
timer.hpp \
trajectory.cpp \
trajectory.hpp \
trajectorydiagnostics.cpp \
trajectorydiagnostics.hpp \
vec3d.hpp \
vectorfield.cpp \
vectorfield.hpp \
verbose.cpp \
verbose.hpp \
gtkfielddiagdialog.cpp \
gtkfielddiagdialog.hpp \
gtkfielddiagwindow.cpp \
gtkfielddiagwindow.hpp \
gtkparticlediagdialog.cpp \
gtkparticlediagdialog.hpp \
gtkparticlediagwindow.cpp \
gtkparticlediagwindow.hpp \
xygraph.cpp \
xygraph.hpp

libibsimu_la_includedir=$(includedir)/ibsimu
libibsimu_la_include_HEADERS = \
fielddiagplot.hpp \
fielddiagplotter.hpp \
particlediagplotter.hpp \
particlediagplot.hpp \
geomplot.hpp \
geomplotter.hpp \
plotter.hpp \
geomplot.hpp \
statusprint.hpp \
fonts.hpp \
fieldgraph.hpp \
label.hpp \
types.hpp \
interpolation.hpp \
gtkplotter.hpp \
ruler.hpp \
histogram.hpp \
legend.hpp \
palette.hpp \
colormap.hpp \
bicgstab.hpp \
bicgstab_solver.hpp \
color.hpp \
ccolmatrix.hpp \
coordmatrix.hpp \
crowmatrix.hpp \
csgobject_solid.hpp \
diag_precond.hpp \
efield.hpp \
empty_precond.hpp \
epot_efield.hpp \
epot_problem.hpp \
error.hpp \
eqpotgraph.hpp \
file.hpp \
frame.hpp \
graph.hpp \
graph3d.hpp \
coordmapper.hpp \
func_solid.hpp \
gauss_seidel.hpp \
geometry.hpp \
hbio.hpp \
ilu0_precond.hpp \
lineclip.hpp \
matrix.hpp \
mvector.hpp \
particledatabase.hpp \
particleiterator.hpp \
particles.hpp \
particlegraph.hpp \
polysolver.hpp \
precond.hpp \
problem.hpp \
qrandom.hpp \
scalarfield.hpp \
scharge.hpp \
scheduler.hpp \
solid.hpp \
solidgraph.hpp \
solver.hpp \
sort.hpp \
timer.hpp \
trajectory.hpp \
trajectorydiagnostics.hpp \
vec3d.hpp \
vectorfield.hpp \
verbose.hpp \
xygraph.hpp

if UMFPACK
libibsimu_la_SOURCES += umfpack_solver.cpp umfpack_solver.hpp
libibsimu_la_include_HEADERS += umfpack_solver.hpp
endif


CLEANFILES = *~


