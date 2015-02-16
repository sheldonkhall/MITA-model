# single-tine-sym-rfa
#
# test for RFA module
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

problemname = "single-tine-sym-rfa"
set_log_level(ERROR) # remove warnings for tests

EM_parameters.cond = 0.132
EM_parameters.V0 = 0.
EM_parameters.Vprobe = 80.
EM_parameters.cond_model = 'nonlinear'
EM_parameters.restrict_mesh = 1
EM_parameters.cond_rate = EM_parameters.cond*0.02
EM_parameters.cond_vap = EM_parameters.cond*0.01

# boundary conditions for electrical problem
EM_parameters.zero_V = [2]
EM_parameters.insulating = [3]
EM_parameters.active_tip = [4]
EM_parameters.symmetry = [1]

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)

# set thermal parameters
thermal_parameters.rho_c_t = 1080.*3455.
thermal_parameters.rho_c_v = 370.*2156. # values for vapourised tissue
thermal_parameters.Lh = 2260.e3 # latent heat of vapourisation
thermal_parameters.Cliq = 0.8 # water content tissue (%)
thermal_parameters.Tu = 373. # upper transition temp
thermal_parameters.Tl = 363. # lower transition temp
thermal_parameters.k = Constant(0.512)
thermal_parameters.dk = Constant(0.02*0.512)
thermal_parameters.omega = Constant(9.19)
thermal_parameters.rho = Constant(1.)
thermal_parameters.c = Constant(3400.)
thermal_parameters.T0 = Constant(310.)
thermal_parameters.T_initial = Constant(310.)
thermal_parameters.restrict_th_mesh = 1  # region to compute thermal solution in
thermal_parameters.bulk_tissue = [2] # fixed temperature (bulk) boundary condition

# solver options
dt_min = .1 # absolute step size minimum (0.0001 good)
dt_max = 1. # absolute step size maximum
tmax = 600 # maximum time (s)
t_out = np.linspace(1,tmax,600) # numpy vector of times at which to save to disk
dt = .1 # time step (s)
thermal_parameters.em_method = 'iterateRFA'
thermal_parameters.k_model = 'linear_limited'
thermal_parameters.perf_model = 'stop'
print "need to check version with nonlinear perfusion, may not be updating"
thermal_parameters.stop_on_me = False

set_log_active(False) # switch off fenics messages

T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)

print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')

