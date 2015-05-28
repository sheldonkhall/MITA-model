# Define domain-specific problemname parameters
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

set_log_level(ERROR) # remove warnings for tests
problemname = "annular-slot-antenna"

EM_parameters.freq = 2.45e9 # Probe frequency
EM_parameters.om = 2 * pi * EM_parameters.freq # Probe angular frequency

# define Power in, inner and outer radius of coaxial dielectric for WPBC
av = 0.01 # applied voltage
EM_parameters.Pin = 30. # power in

# figure 2
#EM_parameters.r_1 = 0.0389
#EM_parameters.r_2 = 0.0974

# figure 6a
EM_parameters.r_1 = 0.00195
EM_parameters.r_2 = 0.00292

# (eps_r, sig, mu_r)
tissue = [[0.0 for x in range(3)] for y in range(3)]
# values from NEVEL paper
tissue[0] = [1.0, 0.0, 1] # air
tissue[1] = [4.51, 0.843*EM_parameters.om*eps_0, 1] # fat / bone
tissue[2] = [49.61, 16.52*EM_parameters.om*eps_0, 1] # muscle

# stupid

# figure 2
# m = 0

# figure 6a
m = 1

EM_parameters.eps_r_by_subdomain = [-1, tissue[0][0], tissue[m][0]]
EM_parameters.sigma_by_subdomain = [-1, tissue[0][1], tissue[m][1]]
EM_parameters.mu_r_by_subdomain = [-1, tissue[0][2], tissue[m][2]]

# WPBC incident field expression from Gas paper
k_dielectric = EM_parameters.om * sqrt(tissue[0][0] * tissue[0][2] * mu_0 * eps_0)
EM_parameters.Z_dielectric = sqrt(mu_0 * tissue[0][2] / (eps_0 * tissue[0][0]))
EM_parameters.C_dielectric = av/ln(r_2 / r_1)
EM_parameters.H_phi_0_re = Expression("(C / (Z * x[0]))",
                        Z=EM_parameters.Z_dielectric, C=EM_parameters.C_dielectric)
EM_parameters.H_phi_0_im = Constant("0.0")


# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)

# compute SAR
thermal_parameters.restrict_th_mesh = 200
thermal_parameters.T0 = Constant(0.)
U, Q, E_r, E_z = compute_SAR_nl(problemname, mesh, interior, boundaries, EM_parameters, thermal_parameters.T0, thermal_parameters)

file_SAR=File("%s/SAR.pvd" % problemname)
file_SAR << Q

print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')
