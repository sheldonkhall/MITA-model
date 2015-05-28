# transient temperature benchmark Vyas and Rustgi 1992
# an SAR is specified mimicking a laser
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
problemname = "vyas-1992"
set_log_level(ERROR) # remove warnings for tests

thermal_parameters.restrict_th_mesh = 1  # region to compute thermal solution in

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)
# dump meshfile
File("mesh/%s.pvd" % problemname) << interior

# set thermal parameters
thermal_parameters.rho_c_t = .0038*1000*1000*1000
thermal_parameters.rho_c_v = .0038*1000*1000*1000
thermal_parameters.Lh = 0. # latent heat of vapourisation
thermal_parameters.Cliq = 0. # water content tissue (%)
thermal_parameters.Tu = 374. # upper transition temp
thermal_parameters.Tl = 372. # lower transition temp
thermal_parameters.Q_sink = 0 # line heat sink
thermal_parameters.k = Constant(4.5837e-04*1000)
thermal_parameters.omega = Constant(0.0013*0.0038*1000*1000*1000)
thermal_parameters.rho = Constant(1.)
thermal_parameters.c = Constant(1.)
thermal_parameters.T0 = Constant(310.)
thermal_parameters.qmet = Constant(0.)
thermal_parameters.rho_t = Constant(1050.)
thermal_parameters.c_t = Constant(3400.)
thermal_parameters.T_initial = Constant(310.)

# solver options
#dt_min = 0.0001 # absolute step size minimum
dt_min = .001
dt_max = .01 # absolute step size maximum
t_out = np.linspace(1,5,5) # numpy vector of times at which to save to disk
dt = .001 # time step (s)
tmax = 5 # maximum time (s)
thermal_parameters.em_method = 'vyas'
thermal_parameters.k_method = 'constant'

T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)

# compute reference solutions
class ref_temp(Expression):
    def eval(self,values,x):

        # map mm to m
        r = x[0]*1000
        z = x[1]*1000
        # set time
        t = 5.

        import scipy.integrate as spint
        import scipy.special as spec
        t0 = 0.3
        E0 = 70*t0
        D = 0.12
        b = 0.0013
        K = 3.5e-4*1e6
        alpha = 0.06
        rhoC = alpha*E0/(pi*t0*K)
        a = 4.5
        tau0 = min([t,t0])
        
        def I(t_):
            y=exp(-b*(t-t_))/(a**2+8*D*(t-t_))*\
               exp(-2*r**2/(a**2+8*D*(t-t_)))*\
               exp(-alpha*z+alpha**2*D*(t-t_))*\
               spec.erfc((2*D*alpha*(t-t_)-z)/(sqrt(4*D*(t-t_))))
            
            return y
            
        y = spint.quad(I,0.,tau0)
        values[0] = K*y[0]+310.
        
comp_ref_temp = ref_temp()
RT = interpolate(comp_ref_temp,T.function_space())

File("%s/analytic_solution.pvd" % problemname) << RT

# def ref_temp(r,z,t):
#     import scipy.integrate as spint
#     import scipy.special as spec
#     t0 = 0.3
#     E0 = 70*t0
#     D = 0.12
#     b = 0.0013
#     K = 3.5e-4*1e6
#     alpha = 0.06
#     rhoC = alpha*E0/(pi*t0*K)
#     a = 4.5
#     tau0 = min([t,t0])

#     def I(t_):
#         y=exp(-b*(t-t_))/(a**2+8*D*(t-t_))*\
#         exp(-2*r**2/(a**2+8*D*(t-t_)))*\
#         exp(-alpha*z+alpha**2*D*(t-t_))*\
#         spec.erfc((2*D*alpha*(t-t_)-z)/(sqrt(4*D*(t-t_))))
        
#         return y

#     y = spint.quad(I,0.,tau0)
#     yk = K*y[0]

#     return yk

print 'start time: ', start_time
print 'end time:   ', tm.strftime('%H:%M:%S')
