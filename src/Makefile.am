lib_LTLIBRARIES = libibsimu.la

libibsimu_la_LDFLAGS = -version-info 0:1:0
#AM_LDFLAGS = -Wall -O3 -march=athlon64
#AM_CPPFLAGS = -Wall -O3 -march=athlon64

AM_LDFLAGS = -Wall -Werror
AM_CPPFLAGS = -Wall -Werror

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
config.hpp \
coordmapper.hpp \
coordmatrix.cpp \
coordmatrix.hpp \
crowmatrix.cpp \
crowmatrix.hpp \
diag_precond.cpp \
diag_precond.hpp \
efield.hpp \
empty_precond.hpp \
epot_efield.cpp \
epot_efield.hpp \
epot_problem.cpp \
epot_problem.hpp \
eqpotgraph.cpp \
eqpotgraph.hpp \
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
matrix.cpp \
matrix.hpp \
meshgraph.cpp \
meshgraph.hpp \
mvector.cpp \
mvector.hpp \
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
timer.cpp \
timer.hpp \
trajectory.cpp \
trajectory.hpp \
trajectorydiagnostics.cpp \
trajectorydiagnostics.hpp \
types.cpp \
types.hpp \
vec3d.hpp \
vectorfield.cpp \
vectorfield.hpp \
verbose.cpp \
verbose.hpp \
xygraph.cpp \
xygraph.hpp

libibsimu_la_includedir=$(includedir)/ibsimu
libibsimu_la_include_HEADERS = \
bicgstab.hpp \
bicgstab_solver.hpp \
ccolmatrix.hpp \
color.hpp \
colormap.hpp \
coordmapper.hpp \
coordmatrix.hpp \
crowmatrix.hpp \
diag_precond.hpp \
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
geomplot.hpp \
geomplotter.hpp \
graph.hpp \
graph3d.hpp \
gs_solver.hpp \
gtkfielddiagexportdialog.hpp \
gtkplotter.hpp \
hbio.hpp \
histogram.hpp \
ilu0_precond.hpp \
interpolation.hpp \
label.hpp \
legend.hpp \
lineclip.hpp \
matrix.hpp \
meshgraph.hpp \
mvector.hpp \
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
types.hpp \
vec3d.hpp \
vectorfield.hpp \
verbose.hpp \
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


