# single-tine-sym-rfa-conv
#
# convergence tests for RFA module
#
# In order to ensure accurate solutions are obtained from the RFA module
# these convergence tests are performed to ensure that the solution is
# independent of the mesh. The solutions that need to be evaluated are
# the electric potential V and the temperature T. The SAR is sensitive
# to the potential, but is essentially the accuracy is limited by this.
# Once the potential has converged adequately, it will be assumed that
# SAR will not improve considerably and will be accepted.
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
import os

start_time = tm.strftime('%H:%M:%S')

problemname = "single-tine-sym-rfa-conv"
set_log_level(ERROR) # remove warnings for tests

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)

# electrical parameters
EM_parameters.cond = 0.2
EM_parameters.V0 = 0.
EM_parameters.Vprobe = 80.
EM_parameters.cond_model = 'nonlinear'
EM_parameters.restrict_mesh = 1
EM_parameters.cond_rate = EM_parameters.cond*0.015
EM_parameters.cond_vap = EM_parameters.cond*0.01

# boundary conditions for electrical problem
EM_parameters.zero_V = [2]
EM_parameters.insulating = [3]
EM_parameters.active_tip = [4]
EM_parameters.symmetry = [1]

# set thermal parameters
# thermal_parameters.rho_c_t = 4.0e6 # normal tissue
# thermal_parameters.rho_c_v = 0.6e6 # values for vapourised tissue
thermal_parameters.rho_c_t = 0. # steady state
thermal_parameters.rho_c_v = 0. # steady state
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
thermal_parameters.restrict_th_mesh = 1  # region to compute thermal solution in
thermal_parameters.bulk_tissue = [2] # fixed temperature (bulk) boundary condition
thermal_parameters.stop_on_me = False
thermal_parameters.cda_update = False

# solver options
dt_min = 1. # absolute step size minimum (0.0001 good)
dt_max = 1. # absolute step size maximum
tmax = 1. # maximum time (s)
t_out = N.array([1.]) # numpy vector of times at which to save to disk
dt = 1. # time step (s)
thermal_parameters.k_model = 'constant'
thermal_parameters.perf_model = 'constant'
thermal_parameters.em_method = 'RFA-const'
thermal_parameters.stop_on_me = False

# look at convergence of a fixed domain
n = 6 # number of times to refine mesh
Qconv=File("%s/Qconv.pvd" % problemname)
Tconv=File("%s/Tconv.pvd" % problemname)

am = np.zeros(n+1)
cm = np.zeros(n+1)
cpow = np.zeros(n+1)
cmin = np.zeros(n+1)
cmax = np.zeros(n+1)
ctmax = np.zeros(n+1)
ctmin = np.zeros(n+1)

for i in range(0,n+1):

    am[i] = str(0.005*0.5**(i))
    cm[i] = str(0.1*0.5**(i))

    # compute SAR
    Q, resistance, power, Vfield = RFA_SAR(problemname, mesh, interior, boundaries, EM_parameters, thermal_parameters.T0, thermal_parameters)
    node_values = Q.vector().array()
    cpow[i] = power
    cmax[i] = node_values.max()
    cmin[i] = node_values.min()
    Qconv << Q

    # compute Temp
    T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)
    Tconv << T
    
    node_values = T.vector().array()
    ctmax[i] = node_values.max()
    ctmin[i] = node_values.min()

    # refine mesh
    meshfile = open('mesh/single-tine-sym-rfa-conv.geo')
    meshfileref = open('mesh/%s-ref.geo' % problemname,'w')
    for line in meshfile:
        # meshfileref.write(line.strip())
        temp = line.replace('am = 0.005','am = %s' % str(0.005*0.5**(i+1)))
        meshfileref.write(temp.replace('cm = 0.1','cm = %s' % str(0.1*0.5**(i+1))))
    meshfile.close()
    meshfileref.close()
    os.system('gmsh -2 mesh/%s-ref.geo -v 0' % problemname)
    os.system('dolfin-convert mesh/%s-ref.msh mesh/%s-ref.xml' % (problemname, problemname))

    # Load new geometry
    mesh = Mesh("mesh/%s-ref.xml" % problemname)
    boundaries = MeshFunction("size_t", mesh, "mesh/%s-ref_facet_region.xml" % problemname)
    interior = MeshFunction("size_t", mesh, "mesh/%s-ref_physical_region.xml" % problemname)
    
print 'antenna mesh: ', am
print 'coarse mesh: ', cm
print 'power: ', cpow
print 'max: ', cmax
print 'min: ', cmin
print 'tmax: ', ctmax
print 'tmin: ', ctmin

# check for domain size convergence
n = 6
for i in range(0,n+1):

    # increase domain
    meshfile = open('mesh/single-tine-sym-rfa-conv.geo')
    meshfileref = open('mesh/%s-ref.geo' % problemname,'w')
    for line in meshfile:
        temp = line.replace('am = 0.005','am = 0.0006')
        temp = temp.replace('cm = 0.1','cm = 0.0125')
        meshfileref.write(temp.replace('rd = 0.075','rd = %s' % str(0.04+0.01*i)))
    meshfile.close()
    meshfileref.close()
    os.system('gmsh -2 mesh/%s-ref.geo -v 0' % problemname)
    os.system('dolfin-convert mesh/%s-ref.msh mesh/%s-ref.xml' % (problemname, problemname))

    # Load new geometry
    mesh = Mesh("mesh/%s-ref.xml" % problemname)
    boundaries = MeshFunction("size_t", mesh, "mesh/%s-ref_facet_region.xml" % problemname)
    interior = MeshFunction("size_t", mesh, "mesh/%s-ref_physical_region.xml" % problemname)

    # compute SAR
    print "running radius: ", str(0.04+0.01*i)
    Q, resistance, power, Vfield = RFA_SAR(problemname, mesh, interior, boundaries, EM_parameters, thermal_parameters.T0, thermal_parameters)
    node_values = Q.vector().array()
    cpow[i] = power
    cmax[i] = node_values.max()
    cmin[i] = node_values.min()
    Qconv << Q

    # compute Temp
    T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)
    Tconv << T
    
    node_values = T.vector().array()
    ctmax[i] = node_values.max()
    ctmin[i] = node_values.min()

print 'power min - max: ', cmin, cmax
print 'temp min - max: ', ctmin, ctmax

print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')

