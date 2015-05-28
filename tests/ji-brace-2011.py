# geometry specified in Ji and Brace 2011
#
# This is really a qualitative verification and validation against another code
# and data. Test the sigmoidal temperature dependence of dielectric properties
# and standard Pennes equation.
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
# define problem name
problemname = "ji-brace-2011"
set_log_level(ERROR) # remove warnings for tests

EM_parameters.freq = 2.45e9 # Probe frequency
EM_parameters.om = 2 * pi * EM_parameters.freq # Probe angular frequency

# define Power in, inner and outer radius of coaxial dielectric for WPBC
EM_parameters.Pin = 25. #######   change all to power and convert here

EM_parameters.r_1 = 0.000255
EM_parameters.r_2 = 0.00084

# materials
# (eps_r, sig, mu_r)
tissue = [[0.0 for x in range(3)] for y in range(4)]

tissue[0] = [1,1,1] # liver dummy
tissue[1] = [2.05,0,1] # di-electric

# stupid way to rearrange
n = 0
m = 1

EM_parameters.eps_r_by_subdomain = [-1, tissue[n][0], tissue[m][0]]
EM_parameters.sigma_by_subdomain = [-1, tissue[n][1], tissue[m][1]]
EM_parameters.mu_r_by_subdomain = [-1, tissue[n][2], tissue[m][2]]

thermal_parameters.restrict_th_mesh = 1  # region to compute thermal solution in

# WPBC incident field expression from Gas paper
# k_dielectric = om * sqrt(tissue[0][0] * tissue[0][2] * mu_0 * eps_0) # unused
EM_parameters.Z_dielectric = sqrt(mu_0 * tissue[1][2] / (eps_0 * tissue[1][0]))
EM_parameters.C_dielectric = sqrt((EM_parameters.Z_dielectric * EM_parameters.Pin) / (pi * ln(EM_parameters.r_2 / EM_parameters.r_1)))
EM_parameters.H_phi_0_re = Expression("(C / (Z * x[0]))",
                        Z=EM_parameters.Z_dielectric, C=EM_parameters.C_dielectric)
EM_parameters.H_phi_0_im = Constant("0.0")


# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)
# dump meshfile
File("mesh/%s.pvd" % problemname) << interior

# values from Ji and Brace
# set thermal parameters
thermal_parameters.rho_c_t = 1050.*3400.
thermal_parameters.rho_c_v = 1050.*3400.
thermal_parameters.Lh = 0. # latent heat of vapourisation
thermal_parameters.Cliq = 0. # water content tissue (%)
thermal_parameters.Tu = 374. # upper transition temp
thermal_parameters.Tl = 372. # lower transition temp
thermal_parameters.Q_sink = 0 # line heat sink
thermal_parameters.k = Constant(0.564)
thermal_parameters.omega = Constant(0.)
thermal_parameters.rho = Constant(1050.)
thermal_parameters.c = Constant(3400.)
thermal_parameters.T0 = Constant(293.)
thermal_parameters.qmet = Constant(0.)
thermal_parameters.rho_t = Constant(1050.)
thermal_parameters.c_t = Constant(3400.)
thermal_parameters.T_initial = Constant(293.)

# set dielectric params
EM_parameters.es1=48.391
EM_parameters.es2=6.286
EM_parameters.es3=0.0764
EM_parameters.es4=1.
EM_parameters.ss1=2.173
EM_parameters.ss2=5.951
EM_parameters.ss3=0.0697
EM_parameters.ss4=0.

# EM_parameters.es1=45.6*2.
# EM_parameters.es2=0.
# EM_parameters.es3=0.
# EM_parameters.es4=0.
# EM_parameters.ss1=1.97*2.
# EM_parameters.ss2=0.
# EM_parameters.ss3=0.
# EM_parameters.ss4=0.

# solver options
dt_min = 0.0001 # absolute step size minimum
#dt_min = 1.
dt_max = 5. # absolute step size maximum
t_out = np.linspace(5,300,60) # numpy vector of times at which to save to disk
dt = 1. # time step (s)
tmax = 300 # maximum time (s)
thermal_parameters.em_method = 'constant'
thermal_parameters.k_model = 'constant'

T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)

print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')
