# Two concentric cylinders of biological tissue with an imposed
# Ez field at fixed radius
# taken from various Paulsen sources
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
problemname = "conc-cyl-test"
set_log_level(ERROR) # remove warnings for tests

EM_parameters.freq = 7.0e7 # Probe frequency
EM_parameters.om = 2 * pi * EM_parameters.freq # Probe angular frequency

# (eps_r, sig, mu_r)
tissue = [[0.0 for x in range(3)] for y in range(12)]
tissue[0] = [10.0, 0.21, 1] # fat
tissue[1] = [70.0, 0.0, 1] # water
tissue[2] = [40.0, 0.35, 1] # lung
tissue[3] = [85.0, 0.802, 1] # muscle
tissue[4] = [1.0, 0.0, 1] # air
tissue[5] = [50.0, 1.44, 1] # bladder
tissue[6] = [10.0, 0.2, 1] # bone
tissue[7] = [113.0, 0.6, 1] # feces
tissue[8] = [80.0, 0.57, 1] # intestine
tissue[9] = [89.0, 1.0, 1] # kidney
tissue[10] = [85.0, 0.802, 1] # tumour
tissue[11] = [89.0, 0.93, 1] # heart

# stupid
n = 2
m = 11

EM_parameters.eps_r_by_subdomain = [-1, tissue[n][0], tissue[m][0]]
EM_parameters.sigma_by_subdomain = [-1, tissue[n][1], tissue[m][1]]
EM_parameters.mu_r_by_subdomain = [-1, tissue[n][2], tissue[m][2]]

EM_parameters.es1=tissue[n][0]*2.
EM_parameters.es2=0.
EM_parameters.es3=0.
EM_parameters.es4=0.
EM_parameters.ss1=tissue[n][1]*2.
EM_parameters.ss2=0.
EM_parameters.ss3=0.
EM_parameters.ss4=0.

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)
# dump meshfile
File("mesh/%s.pvd" % problemname) << interior

thermal_parameters.em_method = 'constant'

# compute SAR
U, Q, E_r, E_z = compute_SAR_nl(problemname, mesh, interior, boundaries, EM_parameters, thermal_parameters.T0, thermal_parameters)

#
# - compute analytic solution
#

from scipy import special

# Get the r and z components
polar = Q.cell().x
r = polar[0]
z = polar[1]

# construct magnitude(E_z)
mag_E_z = project_axisym(sqrt(E_z[0]**2 + E_z[1]**2),Q.function_space())
File("%s/mag-E_z.pvd" % problemname) << mag_E_z

# load variables

ea1 = N.array(0.0,dtype='complex')
ea1 = EM_parameters.eps_r_by_subdomain[1]*eps_0 + 1j*EM_parameters.sigma_by_subdomain[1]/EM_parameters.om
ea2 = N.array(0.0,dtype='complex')
ea2 = EM_parameters.eps_r_by_subdomain[2]*eps_0 + 1j*EM_parameters.sigma_by_subdomain[2]/EM_parameters.om
k1 = EM_parameters.om*N.sqrt(mu_0*ea1)
k2 = EM_parameters.om*N.sqrt(mu_0*ea2)
a = 0.12
R = 0.25

# construct system of equations
A = N.matrix(N.zeros([3,3],dtype='complex'))
b = N.matrix(N.zeros([3,1],dtype='complex'))

A[0,0] = -special.jv(0,k1*a)
A[0,1] = special.jv(0,k2*a)
A[0,2] = special.yv(0,k2*a)

A[1,0] = -k1*special.jv(1,k1*a)
A[1,1] = k2*special.jv(1,k2*a)
A[1,2] = k2*special.yv(1,k2*a)

A[2,0] = 0.0
A[2,1] = special.jv(0,k2*R)
A[2,2] = special.yv(0,k2*R)

b[0] = 0.0
b[1] = 0.0
b[2] = 1.0

# solve system
soln = N.linalg.solve(A,b)

#
# - comparison with analytic solution
#

# Get the r and z components
polar = Q.function_space().cell().x
r = polar[0]
z = polar[1]
    
# project solution onto computational mesh
# various options to compare against analytic. Try something compiled first
exact1 = Constant(soln[0])*bessel_J(0,r)
exact2 = Constant(soln[1])*bessel_J(0,r) + Constant(soln[2])*bessel_Y(0,r)
#exact_1_E_z = project_axisym(exact1,V0_Re) # issue with complex variables here

# overload function object
class exact_E_z(Expression):
    def eval(self,values,x):
        if x[0] < a:
            values[0] = N.abs(soln[0]*special.jv(0,k1*x[0]))
        else:
            values[0] = N.abs(soln[1]*special.jv(0,k2*x[0]) + soln[2]*special.yv(0,k2*x[0]))
f = exact_E_z()
t = interpolate(f,Q.function_space())
File("%s/analytic.pvd" % problemname) << t

# plot solution
#plot(mag_E_z, title="Magnitude of E_z (Numerical solution)")

# L2 norm error
energy = assemble((mag_E_z - t)*dx)
print 'The L2 norm error is: %g' % (energy)

# plot absolute error
h = project_axisym(abs(t - mag_E_z),t.function_space())
File("%s/error.pvd" % problemname) << h
#plot(h, title="Plot of the Absolute Error")
#interactive()

print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')
