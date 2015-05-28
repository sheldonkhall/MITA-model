# empirical SAR from Ai 2012
#
# test case on hold until more general material specification implemented
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
import time as tm

start_time = tm.strftime('%H:%M:%S')

# define problem name
problemname = "ai-2012"
set_log_level(ERROR) # remove warnings for tests

thermal_parameters.restrict_th_mesh = 3  # region to compute thermal solution in

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)

# initialise temperature
thermal_parameters.rho_c_t = 1020.*3628.
thermal_parameters.rho_c_v = 1020.*3628.
thermal_parameters.L = 0 # latent heat of vapourisation
thermal_parameters.Cliq = 0. # water content tissue (%)
thermal_parameters.Tu = 373. # upper transition temp
thermal_parameters.Tl = 353. # lower transition temp
thermal_parameters.Q_sink = 0 # line heat sink
thermal_parameters.k = Constant(0.465)
thermal_parameters.omega = Constant(0.)
thermal_parameters.rho = Constant(1020.)
thermal_parameters.c = Constant(3628.)
thermal_parameters.T0 = Constant(293.)
thermal_parameters.qmet = Constant(0.)
thermal_parameters.k_method = 'ai'
thermal_parameters.em_method = 'ai' # should just do initial
thermal_parameters.T_initial = Constant(293.)

# solver options
#dt_min = 0.0001 # absolute step size minimum
dt_min = 1.
dt_max = 10. # absolute step size maximum
tmax = 200 # maximum time (s)
t_out = N.linspace(1,tmax,10) # numpy vector of times at which to save to disk
dt = 1. # time step (s)

thermal_parameters.perf_model = 'constant'

T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)


print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')
