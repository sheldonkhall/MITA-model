# single-tine-sym-rfa
#
# test for RFA module taken from Trujillo et al. 2012
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

EM_parameters.cond = 0.132
EM_parameters.V0 = 0.
EM_parameters.Vprobe = 80.
EM_parameters.cond_model = 'nonlinear'
EM_parameters.restrict_mesh = 1
EM_parameters.cond_rate = EM_parameters.cond*0.02
EM_parameters.cond_vap = EM_parameters.cond*0.01

# (eps_r, sig, mu_r)
#tissue = [[0.0 for x in range(3)] for y in range(3)]
# values from tosoratti paper
#tissue[0] = [2.1, 0.0, 1] # PTFE coaxial
#tissue[1] = [2.1, 0.0, 1] # PTFE outer sleeve
#tissue[2] = [69, 14.5*EM_parameters.om*eps_0, 1] # muscle

#EM_parameters.eps_r_by_subdomain = [-1, tissue[0][0], tissue[1][0], tissue[2][0]]
#EM_parameters.sigma_by_subdomain = [-1, tissue[0][1], tissue[1][1], tissue[2][1]]
#EM_parameters.mu_r_by_subdomain = [-1, tissue[0][2], tissue[1][2], tissue[2][2]]

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)

#Q = RFA_SAR(problemname, mesh, interior, boundaries, EM_parameters, 310., thermal_parameters)

# set thermal parameters
thermal_parameters.rho_c_t = 1080.*3455.
thermal_parameters.rho_c_v = 370.*2156. # values for vapourised tissue
thermal_parameters.Lh = 2260.e3 # latent heat of vapourisation
thermal_parameters.Cliq = 0.8 # water content tissue (%)
thermal_parameters.Tu = 373. # upper transition temp
thermal_parameters.Tl = 363. # lower transition temp
thermal_parameters.Q_sink = 0 # line heat sink
thermal_parameters.k = Constant(0.512)
thermal_parameters.dk = Constant(0.02*0.512)
thermal_parameters.omega = Constant(9.19)
thermal_parameters.rho = Constant(1.)
thermal_parameters.c = Constant(3400.)
thermal_parameters.T0 = Constant(310.)
thermal_parameters.qmet = Constant(0.)
thermal_parameters.T_initial = Constant(310.)
thermal_parameters.restrict_th_mesh = 1  # region to compute thermal solution in

# solver options
dt_min = .1 # absolute step size minimum (0.0001 good)
dt_max = .1 # absolute step size maximum
tmax = 300 # maximum time (s)
t_out = np.linspace(1,tmax,tmax) # numpy vector of times at which to save to disk
dt = .1 # time step (s)
thermal_parameters.em_method = 'iterateRFA'
thermal_parameters.k_model = 'linear_limited'
thermal_parameters.perf_model = 'stop'
print "need to check version with nonlinear perfusion, may not be updating"
thermal_parameters.stop_on_me = False # (not recommended) switch off check that phase change temp range not exceeded

set_log_active(False) # switch off fenics messages

T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)

print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')

