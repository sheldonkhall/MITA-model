# This script computes the temperature for a volumetric source q constructed
# via the method of moments. No SAR is computed. This is a test of the
# computation of the temperature in axisymmetric coordinates, the symmetry
# condition at the boundary and the Dirichlet condition at the external
# tissue boundaries. For Newton cooling type boundary conditions another
# test would need to be constructed.
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
problemname = "mms-heat-only-test"
set_log_level(ERROR) # remove warnings for tests

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)

# create function space and extract r
V = FunctionSpace(mesh, "CG", order)
# Get the r and z components
rz = V.cell().x
r = rz[0]
z = rz[1]

# set thermal parameters
thermal_parameters.rho_c_t = 0.
thermal_parameters.rho_c_v = 0.
thermal_parameters.k = Constant(0.56)
thermal_parameters.omega = Constant(0.004)
thermal_parameters.rho = Constant(1.)
thermal_parameters.c = Constant(3640.)
thermal_parameters.restrict_th_mesh = 1  # region to compute thermal solution in
thermal_parameters.T0 = Constant(310.)
thermal_parameters.T_initial = Constant(310.)

# define custom source
R = Constant(1.)
eta = sqrt(thermal_parameters.omega*thermal_parameters.c/thermal_parameters.k)
jf = Constant(2.40483)
Lz = Constant(10.)
q = 60*((r*(2*jf**2*Lz**2 + R**2*(pi**2 + 4*Lz**2*eta**2))*bessel_J(0,(jf*r)/R) \
    - 2*jf*Lz**2*(-2*R*bessel_J(1,(jf*r)/R) + jf*r*bessel_J(2,(jf*r)/R))) \
    *cos((pi*z)/(2.*Lz)))/(4.*Lz**2*r*R**2)
thermal_parameters.Q = project_axisym(q*thermal_parameters.k,V)

# solver options
dt_min = 1. # absolute step size minimum (0.0001 good)
dt_max = 1. # absolute step size maximum
tmax = 1. # maximum time (s)
t_out = N.array([1.]) # numpy vector of times at which to save to disk
dt = 1. # time step (s)
thermal_parameters.k_model = 'constant'
thermal_parameters.perf_model = 'constant'
thermal_parameters.em_method = 'custom'
thermal_parameters.stop_on_me = False
thermal_parameters.cda_update = False
thermal_parameters.bulk_tissue = [5] # set correct dirichlet boundary condition
T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)

# This script evaluates the MMS solution
# plot(T)
# interactive()

theta = 60*bessel_J(0,jf*r/R)*cos(pi*z/2/Lz)+310.
temp = theta
T_ref = project_axisym(temp,T.function_space())
T_error = project_axisym(abs(T-T_ref),V)

# L2 norm error
energy = assemble(abs(T-T_ref)*r*dx)
print 'The L2 norm error is: %g' % (energy)

# absolute error plot
File("%s/T-error.pvd" % problemname) << T_error
File("%s/temperature.pvd" % problemname) << project_axisym(T,V)

print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')
