# Comparison against the analytic solution for a coaxial cable
# The first-order absorbing boundary condition performs badly
# can get very accurate solution using wpbc with H_phi_0 = 0
# but needs modified weak form (which is commented out in
# axisymm_mwa module)
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
problemname = "coaxial-test"
set_log_level(ERROR) # remove warnings for tests

# (eps_r, sig, mu_r)
tissue = [[0.0 for x in range(3)] for y in range(1)]
tissue[0] = [2.03, 0, 1] # lossless coaxial

# define Power in, inner and outer radius of coaxial dielectric for WPBC
EM_parameters.Pin = 3.e0
EM_parameters.r_1 = 0.135e-3
EM_parameters.r_2 = 0.470e-3

# stupid
EM_parameters.eps_r_by_subdomain = [-1, tissue[0][0]]
EM_parameters.sigma_by_subdomain = [-1, tissue[0][1]]
EM_parameters.mu_r_by_subdomain = [-1, tissue[0][2]]

#freq = 7.0e7 # Probe frequency
EM_parameters.freq = 2.45e9
EM_parameters.om = 2 * pi * EM_parameters.freq # Probe angular frequency

# WPBC incident field expression from Gas paper
EM_parameters.k_dielectric = EM_parameters.om * sqrt(tissue[0][0] * tissue[0][2] * mu_0 * eps_0)
EM_parameters.Z_dielectric = sqrt(mu_0 * tissue[0][2] / (eps_0 * tissue[0][0]))
EM_parameters.C_dielectric = sqrt((EM_parameters.Z_dielectric * EM_parameters.Pin) / (pi * ln(r_2 / r_1)))
EM_parameters.H_phi_0_re = Expression("(C / (Z * x[0]))",
                        Z=EM_parameters.Z_dielectric, C=EM_parameters.C_dielectric)
EM_parameters.H_phi_0_im = Constant("0.0")

thermal_parameters.restrict_th_mesh = 100 # region to compute thermal solution in

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)

# compute SAR
U, Q, E_r, E_z = compute_SAR_nl(problemname, mesh, interior, boundaries, EM_parameters, thermal_parameters.T0, thermal_parameters)

File("%s/numerical.pvd" % problemname) << U

# comparison against analytic solution
# directly compute the error in the real and imaginary
# parts of the magnetic field.

# construct expressions for analytic H from COMSOL manual (same as Gas)
H_phi_ana = Expression(("C/(x[0]*Z)*cos(k*x[1])","-C/(x[0]*Z)*sin(k*x[1])"),
                       k=EM_parameters.k_dielectric,Z=EM_parameters.Z_dielectric,C=EM_parameters.C_dielectric,pi=pi)
File("%s/analytic.pvd" % problemname) << project_axisym(H_phi_ana,U.function_space())

# plot solution
#plot(U[0], title="Real component of Hphi (Numerical solution)")
#plot(U[1], title="Imaginary component of Hphi (Numerical solution)")
#plot(project_axisym(H_phi_ana[0],Q.function_space()), title="Real component of Hphi (analytic solution)")
#plot(project_axisym(H_phi_ana[1],Q.function_space()), title="Imaginary component of Hphi (analytic solution)")

#plot(sqrt(U[1]**2+U[0]**2), title="Magnitude Hphi (Numerical solution)")
#plot(project_axisym(sqrt(H_phi_ana[0]**2+H_phi_ana[1]**2),Q.function_space()), title="Magnitude Hphi (analytic solution)")

# L2 norm error
energy = assemble((U[0] - H_phi_ana[0])*dx)
print 'The L2 norm error of real part of Hphi is: %g' % (energy)
energy = assemble((U[1] - H_phi_ana[1])*dx)
print 'The L2 norm error of imaginary part of Hphi is: %g' % (energy)

# plot absolute error
h = project_axisym(abs(U - H_phi_ana),U.function_space())
File("%s/error.pvd" % problemname) << h
#plot(h, title="Plot of the Absolute Error")
#plot(h[0], title="Real component of abs(error)")
#plot(h[1], title="Imaginary component of abs(error)")

print 'poor agreement is due to boundary condition not being accurate enough'

#interactive()


print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')
