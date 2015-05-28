# constant temperature cell death test
#
# Copyright (C) 2014  Sheldon Hall (sheldon.hall@eng.ox.ac.uk)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from axisymm_mwa import *
import numpy as np
import time as tm

start_time = tm.strftime('%H:%M:%S')
problemname = "cell-death-test"
set_log_level(ERROR) # remove warnings for tests

thermal_parameters.restrict_th_mesh = 1  # region to compute thermal solution in

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)

# set thermal parameters
thermal_parameters.rho_c_t = 0.
thermal_parameters.rho_c_v = 0.
thermal_parameters.Lh = 0. # latent heat of vapourisation
thermal_parameters.Cliq = 0. # water content tissue (%)
thermal_parameters.Tu = 374. # upper transition temp
thermal_parameters.Tl = 372. # lower transition temp
thermal_parameters.Q_sink = 0 # line heat sink
thermal_parameters.k = Constant(0.)
thermal_parameters.omega = Constant(1.)
thermal_parameters.rho = Constant(1.)
thermal_parameters.c = Constant(1./28.) # should represent 60 C or 333 K
thermal_parameters.T0 = Constant(310.)
thermal_parameters.T_initial = Constant(338.)
thermal_parameters.qmet = Constant(1.)
thermal_parameters.rho_t = Constant(0.)
thermal_parameters.c_t = Constant(0.)

# solver options
#dt_min = 0.0001 # absolute step size minimum
dt_min = 1.
dt_max = 1. # absolute step size maximum
tmax = 900 # maximum time (s)
t_out = np.linspace(1,tmax,tmax) # numpy vector of times at which to save to disk
dt = 1. # time step (s)
thermal_parameters.em_method = 'none'
thermal_parameters.k_method = 'constant'

set_log_active(False) # switch off fenics messages

T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)

# call matlab ode solver for comparison solution
from subprocess import call
call(["matlab","-nodisplay","-r","cd matlab;cell_death_test;exit"])
sol = N.loadtxt("matlab/cell_death_test_sol.dat",delimiter=",")

# store solution in file for paraview
file_matlab=File("%s/matlab.pvd" % problemname)
for t in range(0,tmax):
    T_ref = interpolate(Constant(sol[t,1]),T.function_space())
    T_ref.rename('cell-death',T_ref.label())
    file_matlab << T_ref

print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')
