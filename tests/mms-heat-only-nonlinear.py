# mms-heat-only-nonlinear
#
# uses MMS to test the solution of nonlinear Pennes equation with T dependent k
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

problemname = "mms-heat-only-nonlinear"
set_log_level(ERROR) # remove warnings for tests

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)

# set thermal parameters
thermal_parameters.rho_c_t = 1.
thermal_parameters.rho_c_v = 1.
thermal_parameters.Tu = 373. # upper transition temp
thermal_parameters.Tl = 363. # lower transition temp
thermal_parameters.k = Constant(0.512)
thermal_parameters.dk = Constant(0.02*0.512)
thermal_parameters.omega = Constant(0.)
thermal_parameters.T0 = Constant(310.)
thermal_parameters.T_initial = Expression("310+23*cos(x[0])*sin(x[1])") # initial temperature profile
thermal_parameters.restrict_th_mesh = 1  # region to compute thermal solution in

# solver options
dt_min = .001 # absolute step size minimum (0.0001 good)
dt_max = .001 # absolute step size maximum
tmax = 2 # maximum time (s)
t_out = np.linspace(.2,tmax,10) # numpy vector of times at which to save to disk
dt = .001 # time step (s)
thermal_parameters.k_model = 'linear'
thermal_parameters.perf_model = 'constant'
thermal_parameters.em_method = 'mms-nonlinear'
thermal_parameters.stop_on_me = False
thermal_parameters.cda_update = False
thermal_parameters.bulk_tissue = [5] # set correct dirichlet boundary condition

#set_log_active(False) # switch off fenics messages

T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)

afile = File("%s/mms.pvd" % problemname)

# create ufl function for analytic solution
M = 23.
P = 1.
L = 1.
F = 1.
X = 310.
T_mms = Expression("M*cos(P*x[0])*sin(L*x[1])*exp(-F*t)+X",M=M,P=P,L=L,F=F,X=X,t=0.)
for i in t_out:
    T_mms.t=i
    T_mms_out=interpolate(T_mms,T.function_space())
    T_mms_out.rename('Temperature',T_mms_out.label())
    afile << (T_mms_out,i)

print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')

# L2 norm error
# Get the r and z components
rz = T.function_space().cell().x
r = rz[0]
z = rz[1]

# L2 norm error
energy = assemble(abs(T-T_mms)*r*dx)
print 'The L2 norm error is: %g' % (energy)

# absolute error plot
T_error = project_axisym(abs(T-T_mms),T.function_space())
File("%s/T-error.pvd" % problemname) << T_error
