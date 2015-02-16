# solidification due to line heat sink, 'Heat Conduction' Hahn & Ozisik, 3rd edition p472
# NB: this problem appears to be numerically very difficult to solve. This is due to the temperature tending to -infinity as r->0. This is theoretically fine and the use of quadrature removes evaluating at r=0, but the gradient is very steep and the domain needs to be very long and thin with extremely small elements and time steps. Works as pure verification of equations, but need a better case for testing validity during MWA. The first few time steps violate the restrictions on dT for the apparent heat capacity method and hence stop_on_me must be set to false, but this does not appear to impact final solution.
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
import scipy as sp
from scipy import special, optimize

problemname = "1d-solidification-test"
set_log_level(ERROR) # remove warnings for tests

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)
# dump meshfile
File("mesh/%s.pvd" % problemname) << interior

# restrict mesh to whole domain
thermal_parameters.restrict_th_mesh = 9

# set thermal parameters
thermal_parameters.rho_c_t = 1000.*1762
thermal_parameters.rho_c_v = 1000.*1762
thermal_parameters.Lh = 338e3 # latent heat of vapourisation
thermal_parameters.Cliq = 1. # water content tissue (%)
thermal_parameters.Tu = 274. # upper transition temp
thermal_parameters.Tl = 272. # lower transition temp
thermal_parameters.Q_sink = -1e4 # line heat sink
thermal_parameters.k = Constant(2.22)
thermal_parameters.omega = Constant(0.)
thermal_parameters.rho = Constant(1000.)
thermal_parameters.c = Constant(1.762e3)
thermal_parameters.T0 = Constant(310.)
thermal_parameters.qmet = Constant(0.)
thermal_parameters.rho_t = Constant(1000.)
thermal_parameters.c_t = Constant(1.762e3)
thermal_parameters.T_initial = Constant(310.)
thermal_parameters.stop_on_me = False

# solver options
dt_min = 0.0001 # absolute step size minimum
dt_max = .2 # absolute step size maximum
t_out = np.linspace(0.2,1.,5) # numpy vector of times at which to save to disk
dt = .0001 # time step (s)
tmax = 1. # maximum time (s)
thermal_parameters.em_method = 'none'
thermal_parameters.k_method = 'constant'

T = compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)

# analytic solution

pi = sp.pi

Ti = 310.
Tm = 273. # phase change temperature
ks = 2.22
alphas = ks/1.762e6 # diffusivity
alphal = alphas
rhos = 1000.

# define the transcendental equation
def trans_fun(lam):
    return -thermal_parameters.Q_sink/4./pi*sp.exp(-lam**2) +\
          ks*(Tm-Ti)/sp.special.expn(1,lam**2*alphas/alphal)*\
          sp.exp(-lam**2*alphas/alphal) - lam**2*alphas*rhos*thermal_parameters.Lh

# find the root
lam = sp.optimize.newton(trans_fun,0) # x0 specific for this case

# interface location vs time
def s(t):
    return 2*lam*sp.sqrt(alphas*t)

# solid temperature r=0 -> ri
def Ts(r,t):
    return Tm + -thermal_parameters.Q_sink/4./pi/ks*(sp.special.expn(1,lam**2)-\
        sp.special.expn(1,r**2/4./alphas/t))

# liquid temperature ri -> r infinity
def Tl(r,t):
    return Ti + (Tm - Ti)/sp.special.expn(1,lam**2*alphas/alphal)*\
        sp.special.expn(1,r**2/4./alphal/t)


# create linear function space
w = FunctionSpace(mesh,'CG',1)
Ta = Function(w)

t0 = 0.
class Ta_comp(Expression):
    def eval(self, values, x):
        s0 = s(t0)
        if x[0] <= s0:
            values[0] = Ts(x[0],t0)
        elif x[0] > s0:
            values[0] = Tl(x[0],t0)
        if np.isinf(values[0]):
            values[0] = thermal_parameters.T0
Ta_f = Ta_comp()

afile = File("%s/analytic.pvd" % problemname)

# create ufl function for analytic solution
for i in t_out:
    t0 = i
    Ta = interpolate(Ta_f,w)
    Ta.rename('Temperature',Ta.label())
    afile << Ta

##### need to fix above
