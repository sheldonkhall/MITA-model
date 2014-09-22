# probe specified in Piotr Gas paper
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

problemname = "gas-probe-sensitivity"

EM_parameters.freq = 2.45e9 # Probe frequency
EM_parameters.om = 2 * pi * EM_parameters.freq # Probe angular frequency
EM_parameters.Pin = 30. # power in
EM_parameters.r_1 = 0.000135 # inner conductor
EM_parameters.r_2 = 0.000470 # outer conductor

# materials
# (eps_r, sig, mu_r)
tissue = [[0.0 for x in range(3)] for y in range(4)]

tissue[0] = [43.3,1.69,1] # liver tissue
tissue[1] = [2.03,0,1] # di-electric
tissue[2] = [2.60,0,1] # catheter
tissue[3] = [1,0,1] # air slot

# stupid way to rearrange
n = 0
m = 1
l = 3
p = 2

EM_parameters.eps_r_by_subdomain = [-1, tissue[n][0], tissue[m][0], tissue[l][0], tissue[p][0]]
EM_parameters.sigma_by_subdomain = [-1, tissue[n][1], tissue[m][1], tissue[l][1], tissue[p][1]]
EM_parameters.mu_r_by_subdomain = [-1, tissue[n][2], tissue[m][2], tissue[l][2], tissue[p][2]]

thermal_parameters.restrict_th_mesh = 1  # region to compute thermal solution in

#EM_parameters.es1=tissue[0][0]*2.
#EM_parameters.es2=0.
#EM_parameters.es3=0.
#EM_parameters.es4=0.
#EM_parameters.ss1=tissue[0][1]*2.
#EM_parameters.ss2=0.
#EM_parameters.ss3=0.
#EM_parameters.ss4=0.

EM_parameters.es1=48.391
EM_parameters.es2=6.286
EM_parameters.es3=0.0764
EM_parameters.es4=1.
EM_parameters.ss1=2.173
EM_parameters.ss2=5.951
EM_parameters.ss3=0.0697
EM_parameters.ss4=0.

# WPBC incident field expression from Gas paper
EM_parameters.Z_dielectric = sqrt(mu_0 * tissue[1][2] / (eps_0 * tissue[1][0]))
EM_parameters.C_dielectric = sqrt((EM_parameters.Z_dielectric * EM_parameters.Pin) / (pi * ln(r_2 / r_1)))
EM_parameters.H_phi_0_re = Expression("(C / (Z * x[0]))",
                        Z=EM_parameters.Z_dielectric, C=EM_parameters.C_dielectric)
EM_parameters.H_phi_0_im = Constant("0.0")

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)
# dump meshfile
File("mesh/%s.pvd" % problemname) << interior

# set thermal parameters
thermal_parameters.rho_c_t = 1050.*3628.
thermal_parameters.rho_c_v = 1050.*3628. # values for vapourised tissue
thermal_parameters.Lh = 2260.e3 # latent heat of vapourisation
thermal_parameters.Cliq = 0.8 # water content tissue (%)
thermal_parameters.Tu = 378. # upper transition temp
thermal_parameters.Tl = 368. # lower transition temp
thermal_parameters.Q_sink = 0 # line heat sink
thermal_parameters.k = Constant(0.512)
thermal_parameters.omega = Constant(9.19)
thermal_parameters.rho = Constant(1.)
thermal_parameters.c = Constant(3400.)
thermal_parameters.T0 = Constant(310.)
thermal_parameters.qmet = Constant(0.)
thermal_parameters.T_initial = Constant(310.)

# solver options
dt_min = 1. # absolute step size minimum
#dt_min = 4.
dt_max = 1. # absolute step size maximum
t_out = np.linspace(1,100,100) # numpy vector of times at which to save to disk
dt = 1. # time step (s)
tmax = 100. # maximum time (s)
thermal_parameters.em_method = 'iterate'
thermal_parameters.k_method = 'constant'
thermal_parameters.perf_model = 'stop' # cell state dependent perfusion
thermal_parameters.stop_on_me = False # (not recommended) switch off check that phase change temp range not exceeded
thermal_parameters.cda_update = False

T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)

print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')
