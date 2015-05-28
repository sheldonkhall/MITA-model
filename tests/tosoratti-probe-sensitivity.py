# probe specified in Tosoratti 2003 used in transient calculation.
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

problemname = "tosoratti-probe-sensitivity"
set_log_level(ERROR) # remove warnings for tests

EM_parameters.freq = 2.45e9 # Probe frequency
EM_parameters.om = 2 * pi * EM_parameters.freq # Probe angular frequency
EM_parameters.Pin = 30. # power in
EM_parameters.r_1 = 0.0001016 # inner conductor
EM_parameters.r_2 = 0.000305 # outer conductor

# materials
# (eps_r, sig, mu_r)
tissue = [[0.0 for x in range(3)] for y in range(5)]

# values from tosoratti paper
tissue[1] = [2.2, 0.0, 1] # dielectric
tissue[2] = [10, 0.0, 1] # alumina
tissue[0] = [43.03, 1.69, 1] # liver
tissue[3] = [53.57, 1.81, 1] # muscle
tissue[4] = [69, 1.97, 1] # egg white

m=2

EM_parameters.eps_r_by_subdomain = [-1, tissue[0][0], tissue[1][0], tissue[m][0]]
EM_parameters.sigma_by_subdomain = [-1, tissue[0][1], tissue[1][1], tissue[m][1]]
EM_parameters.mu_r_by_subdomain = [-1, tissue[0][2], tissue[1][2], tissue[m][2]]

thermal_parameters.restrict_th_mesh = 1  # region to compute thermal solution in

# new
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

# compute SAR
thermal_parameters.T0 = Constant(0.)
U, Q, E_r, E_z = compute_SAR_nl(problemname, mesh, interior, boundaries, EM_parameters, thermal_parameters.T0, thermal_parameters)

file_SAR=File("%s/SAR.pvd" % problemname)
file_SAR << Q

# set thermal parameters
thermal_parameters.rho_c_t = 4.0e6 # normal tissue
thermal_parameters.rho_c_v = 0.6e6 # values for vapourised tissue
thermal_parameters.Lh = 2260.e3 # latent heat of vapourisation
thermal_parameters.Cliq = 0.74 # water content tissue (%)
thermal_parameters.Tu = 373. # upper transition temp
thermal_parameters.Tl = 363. # lower transition temp
thermal_parameters.k = Constant(0.5)
thermal_parameters.dk = Constant(0.0033)
thermal_parameters.omega = Constant(50000)
thermal_parameters.rho = Constant(1.)
thermal_parameters.c = Constant(1.)
thermal_parameters.T0 = Constant(310.)
thermal_parameters.T_initial = Constant(310.)
thermal_parameters.bulk_tissue = [5] # fixed temperature (bulk) boundary condition
thermal_parameters.nu = 100 # number of iterations before updating SAR

# solver options
dt_min = 0.001 # absolute step size minimum
#dt_min = 4.
dt_max = .05 # absolute step size maximum
t_out = np.linspace(0.1,50,500) # numpy vector of times at which to save to disk
dt = 0.01 # time step (s)
tmax = 50 # maximum time (s)
thermal_parameters.em_method = 'iterate'
thermal_parameters.k_method = 'constant'
thermal_parameters.perf_model = 'stop' # cell state dependent perfusion

# CHANGE
thermal_parameters.stop_on_me = False
thermal_parameters.cda_update = False

# transient calc
T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)

print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')

