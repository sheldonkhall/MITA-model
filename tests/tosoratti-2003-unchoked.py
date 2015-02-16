# Define domain-specific problemname parameters
# Tosoratti 2003
# This is the probe taken from the paper:
#     A Coaxial Antenna With Miniaturized Choke for Minimally Invasive Interstitial Heating
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

problemname = "tosoratti-2003-unchoked"
set_log_level(ERROR) # remove warnings for tests

EM_parameters.freq = 2.45e9 # Probe frequency
EM_parameters.om = 2 * pi * EM_parameters.freq # Probe angular frequency

# define Power in, inner and outer radius of coaxial dielectric for WPBC
#av = 1 # applied voltage
EM_parameters.Pin = 50./2/N.pi # applied power

# radius for WPBC
EM_parameters.r_1 = 0.0001435
EM_parameters.r_2 = 0.00047

# (eps_r, sig, mu_r)
tissue = [[0.0 for x in range(3)] for y in range(2)]
# values from tosoratti paper
tissue[0] = [2.1, 0.0, 1] # PTFE coaxial
tissue[1] = [69, 14.5*EM_parameters.om*eps_0, 1] # muscle

EM_parameters.eps_r_by_subdomain = [-1, tissue[1][0], tissue[0][0]]
EM_parameters.sigma_by_subdomain = [-1, tissue[1][1], tissue[0][1]]
EM_parameters.mu_r_by_subdomain = [-1, tissue[1][2], tissue[0][2]]

thermal_parameters.restrict_th_mesh = 1 # region to compute thermal solution in

# WPBC incident field expression from Gas paper
EM_parameters.Z_dielectric = sqrt(mu_0 * tissue[0][2] / (eps_0 * tissue[0][0]))
#EM_parameters.C_dielectric = av/ln(r_2 / r_1)
EM_parameters.C_dielectric = sqrt((EM_parameters.Z_dielectric * EM_parameters.Pin) / (pi * ln(r_2 / r_1)))
EM_parameters.H_phi_0_re = Expression("(C / (Z * x[0]))",
                        Z=EM_parameters.Z_dielectric, C=EM_parameters.C_dielectric)
EM_parameters.H_phi_0_im = Constant("0.0")


# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)

# set thermal parameters
thermal_parameters.rho_c_t = 1060.*3411.
thermal_parameters.rho_c_v = 1060.*3411.
thermal_parameters.Lh = 0. # latent heat of vapourisation
thermal_parameters.Cliq = 0. # water content tissue (%)
thermal_parameters.Tu = 374. # upper transition temp
thermal_parameters.Tl = 372. # lower transition temp
thermal_parameters.Q_sink = 0 # line heat sink
thermal_parameters.k = Constant(0.56)
thermal_parameters.omega = Constant(0.004)
thermal_parameters.rho = Constant(1020.)
thermal_parameters.c = Constant(3640.)
thermal_parameters.T0 = Constant(310.)
thermal_parameters.qmet = Constant(0.)
thermal_parameters.rho_t = Constant(1060.)
thermal_parameters.c_t = Constant(3411.)

# set dielectric params
EM_parameters.es1=tissue[1][0]*2.
EM_parameters.es2=0.
EM_parameters.es3=0.
EM_parameters.es4=0.
EM_parameters.ss1=tissue[1][1]*2.
EM_parameters.ss2=0.
EM_parameters.ss3=0.
EM_parameters.ss4=0.

# solver options
#dt_min = 0.0001 # absolute step size minimum
dt_min = .1
dt_max = 10. # absolute step size maximum
t_out = N.linspace(5,300,60) # numpy vector of times at which to save to disk
dt = .5 # time step (s)
tmax = 300 # maximum time (s)

# compute SAR
U, Q, E_r, E_z = compute_SAR_nl(problemname, mesh, interior, boundaries, EM_parameters, thermal_parameters.T0, thermal_parameters)

file_SAR=File("%s/SAR.pvd" % problemname)
file_SAR << Q
