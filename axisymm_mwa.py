## Python Module: axisymm_mwa
## Author: Sheldon Hall
## Email: sheldon.hall@eng.ox.ac.uk
##
## This module contains functions for predicting the outcomes of
## minimally invasive cancer treatments (MICTs). Namely: RFA,
## MWA, CA, IRE. The main focus of this code is to perform
## sensitivity analyses on representative (simplified) problems.
##
## The bioheat model chosen is the effective heat capacity form
## of Pennes equations, which utilises the only transient
## solver. The computations of SAR for RFA and MWA are quasi-
## static and performed as required by the nonlinear solvers for
## the bioheat equation.
##
## Several classes are included to define data structures for
## various model, solver and code parameters.

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

#
# - includes
#

from scipy import constants as S
import numpy as N
from dolfin import *
import sys
import os

#
# - set fenics params
#

parameters["num_threads"] = 6 # invoke multi-thread parallel support
parameters["allow_extrapolation"] = True # for mesh mappings
parameters["form_compiler"]["optimize"] = True # optimise compiler pg 167 fenics book

#solver = NewtonSolver("mumps")
#solver = LUSolver("petsc")
#solver = KrylovSolver("cg","ilu")

#
# - Some constants
#

pi = S.pi
e = S.e
c0 = S.c
mu_0 = S.mu_0
eps_0 = S.epsilon_0

#
# - dummies to catch errors
#

Pin = 0. # input power for WPBC
r_1 = 1. # inner radius for coaxial feed
r_2 = 2. # outer radius for coaxial feed
H_phi_0_re = Constant(0.0) # input real amplitude for WPBC
H_phi_0_im = Constant(0.0) # input imaginary amplitude for WPBC
qmet = Constant(0.) # metabolic heat term

#
# - global parameters
#

order = 2 # order of elements

#
# - EM parameters
#

class EM_parameters:

    restrict_mesh = 1

    # RFA parameters
    cond = 0.333
    cond_rate = 0.
    V0 = 0. # ground pad voltage
    Vprobe = 100. # probe voltage
    cond_model = 'constant' # dependence of electrical conductivity
    cond_vap = 0.5 # conductivity after vapourisation
    imp_max = 120 # impedance threshold for control system
    imp_t_off = 15 # time to switch off probe after imp_max exceeded
    mms_source = Constant(0.) # source term for use in mms

    # MWA parameters
    freq = 2.45e9 # Probe frequency
    om = 2 * pi * freq # Probe angular frequency
    Pin = 0. # input power for WPBC
    r_1 = 1. # inner radius for coaxial feed
    r_2 = 2. # outer radius for coaxial feed
    H_phi_0_re = Constant(0.0) # input real amplitude for WPBC
    H_phi_0_im = Constant(0.0) # input imaginary amplitude for WPBC
    eps_r_by_subdomain = [] # baseline epsilon r values by subdomain
    sigma_by_subdomain = [] # baseline sigma values by subdomain
    mu_r_by_subdomain = [] # baseline mu r values by subdomain
    Z_dielectric = 0. # WPBC
    C_dielectric = 0. # WPBC

    # temperature dependent dielectric
    es1=48.391
    es2=6.286
    es3=0.0764
    es4=1.
    ss1=2.173
    ss2=5.951
    ss3=0.0697
    ss4=0.

    # boundaries
    zero_V = []
    insulating = []
    symmetry = []
    active_tip = []
    
#    def __init__(self,Pin):
#        self.Pin=Pin

#
# - Thermal parameters
#

class thermal_parameters:
    restrict_th_mesh = 1  # region to restrict thermal solution to
    qmet = Constant(0.) # metabolic heat term
    rho_c_t = 1060.*3411. # rho*c in tissue (phase)
    rho_c_v = 4.4e5 # rho*c in vapourised tissue (phase)
    Lh = 0 # latent heat of vapourisation
    Cliq = 0. # water content tissue (%)
    Tu = 374. # upper transition temp
    Tl = 372. # lower transition temp
    Q_sink = 0 # line heat sink
    k = Constant(0.56) # thermal conductivity
    dk = 0. # rate of change of thermal conductivity
    omega = Constant(0.004) # blood perfusion
    rho = Constant(1020.) # density blood
    c = Constant(3640.) # specific heat blood
    T0 = Constant(310.) # baseline temperature
    k_model = 'constant' # choose model for thermal conductivity
    em_method = 'constant' # choose why type of EM model to use
    T_initial = Constant(310.) # initial flat temperature profile
    perf_model = 'constant'
    stop_on_me = True # error given when T change over dt too large
    Q = Constant(0.) # allow custom heat source to be given
    cda_update = True # compute cell death
    p_stop = 0.8 # value of viability at which to stop perfusion
    cool_probe_temp = 310. # coolant temp on probe boundary
    cool_probe = 100. # boundary to apply condition
    h_transfer = 0. # heat transfer coefficient into coolant
    nu = 1. # number of iterations of heat solver before updating heat source

    # boundaries
    bulk_tissue = []
    cool_probe = []

class cell_death_parameters:
    kb = 7.77e-3
    kf_ = 3.33e-3
    Tk = 40.5
    A_init = 0.99

#
# - Functions
#

# - create output directory
def ensure_dir(f):
    if not os.path.exists(f):
        os.makedirs(f)

# compute_SAR_nl
#
# takes the temperature field, mesh and material properties as arguments and returns SAR
#
# arguments:
#
# returns:

#def compute_SAR_nl(mesh, interior, boundaries, problemname, eps_r_by_subdomain,
#                   mu_r_by_subdomain, sigma_by_subdomain, H_phi_0_re, H_phi_0_im, om,
#                   es1, es2, es3, es4, ss1, ss2, ss3, ss4, T):
def compute_SAR_nl(problemname, mesh, interior, boundaries, emp, T, thp):

    # check directory exists for results
    ensure_dir(problemname)

    # set solver params in advance
    #solver = KrylovSolver("cg", "hypre_euclid")

    # set measure
    dss = ds[boundaries]
    
    #
    # - Define the function spaces
    #

    V0 = VectorFunctionSpace(mesh, "CG", order, dim=2) # complex scalar field
    V0E = VectorFunctionSpace(mesh, "CG", order, dim=4) # complex vector field
    V0_Re = FunctionSpace(mesh, "CG", order) # scalar field
    V0_dc = FunctionSpace(mesh, "DG", 0) # material properties (discontinuous on boundary)

    #
    # - Define piecewise constant material properties
    #

    # vectorised numpy quicker than a python loop
    eps_r = Function(V0_dc)
    eps_r.vector()[:] = N.choose(N.asarray(interior.array(), dtype=N.int32), emp.eps_r_by_subdomain)

    mu_r = Function(V0_dc)
    mu_r.vector()[:] = N.choose(N.asarray(interior.array(), dtype=N.int32), emp.mu_r_by_subdomain)

    sigma = Function(V0_dc)
    sigma.vector()[:] = N.choose(N.asarray(interior.array(), dtype=N.int32), emp.sigma_by_subdomain)

    # substitue values in tissue
    # take restrict_th_mesh and set T dependent properties in that region
    T_p = Function(V0_dc)
#    T_p = interpolate(T,V0_dc)
    T_p = project_axisym(T,V0_dc)
    T_array = T_p.vector().array()
    eps_r_array = eps_r.vector().array()
    eps_r_array[interior.array()==thp.restrict_th_mesh] = emp.es1*(1-1/(1+N.exp(emp.es2-emp.es3*(T_array[interior.array()==thp.restrict_th_mesh]+37.))))+emp.es4
    eps_r.vector()[:] = eps_r_array
    sigma_array = sigma.vector().array()
    sigma_array[interior.array()==thp.restrict_th_mesh] = emp.ss1*(1-1/(1+N.exp(emp.ss2-emp.ss3*(T_array[interior.array()==thp.restrict_th_mesh]+37.))))+emp.ss4
    sigma.vector()[:] = sigma_array
    #eps_r = es1*(1-1/(1+exp(es2-es3*(T-273.))))+es4
    #sigma = ss1*(1-1/(1+exp(ss2-ss3*(T-273.))))+ss4

    File("%s/eps_r.pvd" % problemname) << eps_r
    File("%s/sigma.pvd" % problemname) << sigma

    #
    # - construct weak form and boundary conditions
    #

    # Finite element test and trial
    H_phi = TrialFunction(V0)
    T = TestFunction(V0)

    # Get the r and z components
    polar = V0.cell().x
    r = polar[0]
    z = polar[1]

    # Get the surface normal
    n = V0.cell().n

    ### should be moved to specific example
    # imposed electric field for BC
    E0 = Constant(1.0)
    f = Constant(0.0)

    # define k_0
    k_0 = emp.om * sqrt(eps_0 * mu_0)

    # Reciprocal of complex relative permittivity
    mod_Eps_r = eps_r * eps_r + sigma * sigma / (emp.om * emp.om * eps_0 * eps_0)
    reEps_r_Re = eps_r / mod_Eps_r
    reEps_r_Im = sigma / (emp.om * eps_0) / mod_Eps_r

    # Complex relative permittivity and square root all materials
    # equivalent to mod_Eps_r
    mer = N.array(emp.eps_r_by_subdomain)**2 + (N.array(emp.sigma_by_subdomain)/(emp.om*eps_0))**2
    # equivalent to reEps_r_Re + j*reEps_r_Im above
    rer = N.array(emp.eps_r_by_subdomain)/mer + N.array(emp.sigma_by_subdomain)/(emp.om*eps_0)/mer*1j
    # square root of complex number
    srer = N.sqrt(rer)
    # extract real and imaginary part and assign to dolfin variable
    srer_re = Function(V0_dc)
    srer_re.vector()[:] = N.choose(N.asarray(interior.array(),dtype=N.int32),N.real(srer))
    srer_im = Function(V0_dc)
    srer_im.vector()[:] = N.choose(N.asarray(interior.array(),dtype=N.int32),N.imag(srer))

    #
    # - construct weak form + BCs
    #

    ##     boundaries defined as:
    ##         1 - symmetry condition at x = 0
    ##         2 - first order absorbing boundary condition
    ##         3 - imposed z-component of electric field
    ##         4 - waveguide port boundary condition
    
    # Main operators (applied)
    curl_H_phi_r = - H_phi.dx(1)
    curl_H_phi_z = (1 / r) * (r * H_phi).dx(0)
    curl_T_r = - T.dx(1)
    curl_T_z = (1 / r) * (r * T).dx(0)

    # Define the bilinear forms
    #    Mass form
    m = mu_r * (k_0 ** 2) * inner(T, H_phi)
    s = reEps_r_Re * (curl_T_r[0] * curl_H_phi_r[0] + curl_T_r[1] * curl_H_phi_r[1]) \
        + reEps_r_Im * (curl_T_r[1] * curl_H_phi_r[0] - curl_T_r[0] * curl_H_phi_r[1]) \
        + reEps_r_Re * (curl_T_z[0] * curl_H_phi_z[0] + curl_T_z[1] * curl_H_phi_z[1]) \
        + reEps_r_Im * (curl_T_z[1] * curl_H_phi_z[0] - curl_T_z[0] * curl_H_phi_z[1])

    a = r * s * dx - \
        r * m * dx + \
        -r * k_0 * (T[1] * (reEps_r_Im*H_phi[1]-reEps_r_Re*H_phi[0]) +\
                   T[0] * (reEps_r_Im*H_phi[0] + reEps_r_Re*H_phi[1])) * dss(2) \
        - r * k_0 * T[0] * (srer_im * (H_phi[0]) + srer_re * (H_phi[1])) * dss(4) \
        - r * k_0 * T[1] * (srer_im * (H_phi[1]) - srer_re * (H_phi[0])) * dss(4)

    ### gives better solution to coax problem wpbc without source
    # a = r * s * dx - \
    #     r * m * dx + \
    #     - r * k_0 * T[0] * (srer_im * (H_phi[0]) + srer_re * (H_phi[1])) * dss(2) \
    #     - r * k_0 * T[1] * (srer_im * (H_phi[1]) - srer_re * (H_phi[0])) * dss(2) \
    #     - r * k_0 * T[0] * (srer_im * (H_phi[0]) + srer_re * (H_phi[1])) * dss(4) \
    #     - r * k_0 * T[1] * (srer_im * (H_phi[1]) - srer_re * (H_phi[0])) * dss(4)

    L = r * (T[0]+T[1]) * f * dx + r * emp.om * eps_0 * E0 * T[1] * dss(3) +\
        r * k_0 * T[0] * (srer_im * (- 2*emp.H_phi_0_re) + srer_re * (- 2*emp.H_phi_0_im)) * dss(4) +\
        r * k_0 * T[1] * (srer_im * (- 2*emp.H_phi_0_im) - srer_re * (- 2*emp.H_phi_0_re)) * dss(4)

    
    bc1 = DirichletBC(V0, Constant((0.0, 0.0)), boundaries, 1)
    bcs=[bc1]

    #
    # - solve for H_phi in axisymmetric case
    #

    U = Function(V0)
    #solve(a == L, U, bcs,
    #      solver_parameters={"linear_solver": "mumps",
    #                         "preconditioner": "hypre_euclid"})
    solve(a == L, U, bcs,
          solver_parameters={"linear_solver": "mumps",
                "preconditioner": "hypre_euclid"})

    #
    # - Post-processing
    #

    # compute E_r component
    uE = TrialFunction(V0)
    TE = TestFunction(V0)
    aE = r * inner(TE, uE) * dx
    LE = r * (1 / (emp.om * eps_0)) * \
         ((- reEps_r_Im * U[0].dx(1)- reEps_r_Re * U[1].dx(1)) * TE[0]) * dx + \
         r * (1 / (emp.om * eps_0)) * \
         ((- reEps_r_Im * U[1].dx(1) + reEps_r_Re * U[0].dx(1)) * TE[1]) * dx
    E_r = Function(V0)
    #solve(aE == LE, E_r, solver_parameters={"linear_solver": "mumps"})
    solve(aE == LE, E_r,
          solver_parameters={"linear_solver": "mumps",
                "preconditioner": "hypre_euclid"})

    # compute E_z component
    aE = r * inner(TE, uE) * dx
    LE = r * (1 / (emp.om * eps_0)) * \
         (reEps_r_Im * ((1 / r) * (r * U[0]).dx(0)) + \
         reEps_r_Re * ((1 / r) * (r * U[1]).dx(0))) * TE[0] * dx + \
         r * (1 / (emp.om * eps_0)) * (reEps_r_Im * ((1 / r) * (r * U[1]).dx(0)) - \
         reEps_r_Re * ((1 / r) * (r * U[0]).dx(0))) * TE[1] * dx
    E_z = Function(V0)
    #solve(aE == LE, E_z, solver_parameters={"linear_solver": "mumps"})
    solve(aE == LE, E_z,
          solver_parameters={"linear_solver": "mumps",
                "preconditioner": "hypre_euclid"})

    # compute SAR
    #Q = project_axisym(0.5 * sigma * (E_r[0] ** 2 + E_r[1] ** 2 + E_z[0] ** 2 + E_z[1] ** 2),V0_Re)
    Q = project_axisym(0.5 * sigma * (E_r[0] ** 2 + E_r[1] ** 2 + E_z[0] ** 2 + E_z[1] ** 2),V0_dc)

    # compute power according to RFA
    power = assemble(Q*r*dx)*2*N.pi
    print "power: ", power

    # File("%s/U.pvd" % problemname) << U
    # File("%s/E_r.pvd" % problemname) << E_r
    # File("%s/E_z.pvd" % problemname) << E_z
    # File("%s/Q.pvd" % problemname) << Q

    # normalize solution to 1
    ## U_nodal_values = E_z.vector() # extract nodal values
    ## U_array = U_nodal_values.array() # copy to numpy array
    ## U_max = U_array.max() # numpy find max value
    ## U_array /= U_max
    ## E_z.vector()[:] = U_array
    ## E_z.vector().set_local(U_array) # alternative

    ## By now the solution should have been computed and stored. Anything following
    ## this is problem specific and can just be specified in a separate script to
    ## keep things tidy

    return U, Q, E_r, E_z

# compute_T_enthalpy_nl
#
# solves the time-dependent enthalpy form of the bioheat equation using backward
# differences. This form of the equation allows the computation of phase changes
# in the tissue as a result of heating and freezing. The main extension over the
# most simple time-dependent bioheat equation is the ability to take functions
# of temperature and time as parameters. This therefore also includes nonlinear
# phenomena.
#
# arguments:
#
# returns:

def compute_enthalpy_nl(mesh, interior, boundaries, problemname, dt, tmax, dt_min, dt_max, t_out, thp, emp):

    #checks
    ensure_dir(problemname) # check directory exists for results
    eps = N.finfo(float).eps # useful quantity
    if t_out.size == 1 and t_out[0] > dt_max-eps: # catch output error
        print 'single time point out'
    elif N.any(N.diff(t_out)<dt_max-eps):
        error("largest time step spans more than one output reduce dt_max or coarsen t_out")

    print "--+--+-- start time-dependent bioheat solve --+--+--"        

    # NOTE:
    # 
    # EQUATIONS HAVE BEEN SCALED IN TERMS OF THETA = T - T0
    thp.Tu = thp.Tu - 310
    thp.Tl = thp.Tl - 310

    # set solver params in advance
    solver = KrylovSolver("cg", "hypre_euclid")
    #solver.parameters["absolute_tolerance"] = 1E-7
    #solver.parameters["relative_tolerance"] = 1E-4
    #solver.parameters["maximum_iterations"] = 1000
    #set_log_level(DEBUG)

    # output files
    file_temp=File("%s/enthalpy.pvd" % problemname)
    file_SAR=File("%s/SAR.pvd" % problemname)
    file_cd=File("%s/cell-death.pvd" % problemname)
    file_perf=File("%s/perfusion.pvd" % problemname)
    file_vfield=File("%s/voltage.pvd" % problemname)

    # define a restriction
    # to generalise this define new meshfunction that sets to 1 everything in restriction
    interior_new = MeshFunction("uint",interior)
    help = N.asarray(interior_new.array())
    for ent in N.nditer(help, op_flags=['readwrite']):
        ent[...] = N.where(N.any(ent == thp.restrict_th_mesh),1,0)
    interior_new.array()[:] = help

    restriction = Restriction(interior_new,1) # restrict thermal calc to tissue only
    W = FunctionSpace(restriction, 'CG', order)
    W_dg = FunctionSpace(restriction, 'DG', order) # DG SAR
    W_dg_0 = FunctionSpace(restriction, 'DG', 0) # material properties (perfusion etc)

    # set measure
    dss = ds[boundaries]

    # Get the r and z components
    polar = W.cell().x
    r = polar[0]
    z = polar[1]

    # define quantities that need updating
    dte = Expression('dt',dt=0.)
    cur_time = Expression('t',t=0.)

    # initial uniform temperature
#    T_prev = interpolate(thp.T_initial,W)
    T_prev = project_axisym(thp.T_initial-310,W)

    # initial values (if needed)
    resistance = 0.
    
    print "--+--+-- initial SAR                        --+--+--"        
    # initial SAR
    if thp.em_method=='iterate' or thp.em_method=='constant':
        U, Q, E_r, E_z = compute_SAR_nl(problemname, mesh, interior, boundaries, emp, T_prev, thp)
    elif thp.em_method=='ai':
        Q = interpolate(Constant(1.),W)
    elif thp.em_method=='none':
        Q = interpolate(Constant(0.),W)
    elif thp.em_method=='vyas':
        Q = Expression('2*60*70/pi/pow(0.0045,2)*exp(-2*pow(x[0],2)/pow(0.0045,2)-60*x[1])')
    elif thp.em_method=='RFA-const' or thp.em_method=='iterateRFA':
        Q, resistance, power, Vfield = RFA_SAR(problemname, mesh, interior, boundaries, emp, T_prev, thp)
    elif thp.em_method=='custom':
        Q = thp.Q
    elif thp.em_method=='mms-nonlinear':
        M = 23.
        P = 1.
        L = 1.
        H = 0.512
        R = 0.02*0.512
        F = 1.
        Q = (M*(M*pow(P,2)*r*R*pow(cos(P*r),2) - 2*pow(L,2)*M*r*R*pow(cos(P*r),2)*cos(2*L*z) - M*pow(P,2)*r*R*pow(cos(P*r),2)*cos(2*L*z) + 2*exp(F*cur_time)*H*pow(L,2)*r*cos(P*r)*sin(L*z) + 2*exp(F*cur_time)*H*pow(P,2)*r*cos(P*r)*sin(L*z) - 2*exp(F*cur_time)*F*r*thp.rho_c_t*cos(P*r)*sin(L*z) + 2*exp(F*cur_time)*H*P*sin(P*r)*sin(L*z) - 2*M*pow(P,2)*r*R*pow(sin(P*r),2)*pow(sin(L*z),2) + M*P*R*sin(2*P*r)*pow(sin(L*z),2)))/(2.*exp(2*F*cur_time)*r)
    # elif thp.em_method=='mms-nonlinear-full':
    #     M = 310.
    #     P = 1.
    #     L = 1.
    #     H = 0.512
    #     R = 0.02*0.512
    #     F = 1.
    #     W1 = .15
    #     Y = .02*.15
    #     G = 1.
    #     A1 = 1.
    #     B = 1.
    #     emp.mms_source = Constant(0.)
    #     Q, resistance, power, Vfield = RFA_SAR(problemname, mesh, interior, boundaries, emp, T_prev, thp)
    
    # interpolate heat source onto restriction
#    qext = interpolate(Q,W)
    qext = project_axisym(Q,W)

    # plot(qext)
    # interactive()

    if thp.em_method=='vyas':
        qext = conditional(And(gt(z,0),lt(cur_time,0.3)),Q,0.)
    elif thp.em_method=='mms-nonlinear':
        qext = Q
    # elif thp.em_method=='mms-nonlinear-full':
    #     qext = Q

    # apply boundary conditions according to mesh function
    bcs = []
    for index in thp.bulk_tissue:
        bcs.append(DirichletBC(W, 0., boundaries, index))

    # for index in thp.cool_probe:
    #     bcs.append(DirichletBC(W, thp.cool_probe_temp-310., boundaries, index))

    # for a neumann condition to create a heat sink at r=0 set gmsh to 6

    # define variational problem
    T = TrialFunction(W)
    v = TestFunction(W)
    f = qext
    q = thp.qmet

    # define effective heat capacity using ufl conditional

    rc1 = conditional(lt(T_prev,thp.Tl), thp.rho_c_t, 0.)
    rc2 = conditional(And(ge(T_prev,thp.Tl),le(T_prev,thp.Tu)), (thp.rho_c_t+thp.rho_c_v)/2+thp.rho*thp.Lh*thp.Cliq*(1./(thp.Tu-thp.Tl)), 0.)
    rc3 = conditional(gt(T_prev,thp.Tu), thp.rho_c_v, 0.)

    # define thermal conductivity
    k=thp.k
    if thp.k_model=='linear':
        k=thp.k + thp.dk*(T_prev)
    elif thp.k_model=='linear_limited':
        k=conditional(le(T_prev,63.),thp.k + thp.dk*(T_prev),thp.k + thp.dk*(63.))
    elif thp.k_model=='ai':
        k=conditional(lt(T_prev,43), 0.465, 0.) + conditional(And(ge(T_prev,43),lt(T_prev,73), 0.867, 0.)) + conditional(ge(T_prev,73),1.460,0.)

    # define perfusion term using D > 0.8 as shutoff
    D_prev=interpolate(Constant(0.),W)
    omega=thp.omega
    # project D onto piecewise constant mesh to stop negative values
    D_prev_const = project_axisym(D_prev,W_dg_0)
    if thp.perf_model=='stop':
        omega=conditional(gt(D_prev_const,thp.p_stop), thp.omega, 0.)
        print "check perfusion threshhold"

    # old unscaled and unstable weak form    
    # a = k*inner(nabla_grad(T), nabla_grad(v))*r*dx + v*omega*thp.rho*thp.c*T*r*dx + v*rc1/dte*T*r*dx + v*rc2/dte*T*r*dx + v*rc3/dte*T*r*dx
    # L = f*v*r*dx+q*v*r*dx+v*omega*thp.rho*thp.c*thp.T0*r*dx + v*rc1/dte*T_prev*r*dx + v*rc2/dte*T_prev*r*dx + v*rc3/dte*T_prev*r*dx + v*thp.Q_sink/(2*pi)*dss(6)

    # scaled but no heat transfer coefficient
    # a = k*inner(nabla_grad(T), nabla_grad(v))*r*dx + v*omega*thp.rho*thp.c*T*r*dx + v*rc1/dte*T*r*dx + v*rc2/dte*T*r*dx + v*rc3/dte*T*r*dx
    # L = f*v*r*dx + q*v*r*dx + v*rc1/dte*T_prev*r*dx + v*rc2/dte*T_prev*r*dx + v*rc3/dte*T_prev*r*dx + v*thp.Q_sink/(2*pi)*dss(6)

    # heat transfer
    a = k*inner(nabla_grad(T), nabla_grad(v))*r*dx + v*omega*thp.rho*thp.c*T*r*dx + v*rc1/dte*T*r*dx + v*rc2/dte*T*r*dx + v*rc3/dte*T*r*dx + v*T*thp.h_transfer*r*dss(4)

    L = f*v*r*dx + q*v*r*dx + v*rc1/dte*T_prev*r*dx + v*rc2/dte*T_prev*r*dx + v*rc3/dte*T_prev*r*dx + v*thp.Q_sink/(2*pi)*dss(6)


    T = Function(W)
    T_out = Function(W)

    # set initial temperature for SAR
#    T = interpolate(thp.T_initial,W)
    T = project_axisym(thp.T_initial-310.,W)

    # assemble in advance of time iteration
    A = None
    b = None

    # SAR update criteria
    Q = Function(W_dg)
    iu = thp.nu-1
    store_resistance = N.array(t_out) # save the resistance at output times
    store_power = N.array(t_out) # save power
    power = 0.

    # initialise cell death
    n = len(T.vector().array())
    cda = N.zeros(2*n) # cell death array
    cda[::2] = cell_death_parameters.A_init
    D = interpolate(Constant(0.),W) # dead field
 
    # control system
    imp_on = True
    imp_off_t_start = dt

    t = dt
    step_inc = True
    while t <= tmax+eps:

        # update SAR
        # this is SLOW so introduce index iu
        # iu increments each transient iteration, when it reaches nu it updates SAR
        iu += 1
        if (iu == thp.nu) and (thp.em_method=="iterate"):
            print "--+--+-- updating SAR                       --+--+--"        
            U, Q, E_r, E_z = compute_SAR_nl(problemname, mesh, interior, boundaries, emp, T,thp)
            iu = 0
        if (iu == thp.nu) and (thp.em_method=="iterateRFA"):
            if imp_on: # control system switch
                print "--+--+-- updating SAR                       --+--+--"        
                Q, resistance, power, Vfield = RFA_SAR(problemname, mesh, interior, boundaries, emp, T, thp)
                iu = 0

            else:
                print "--+--+-- SAR offfffffff                     --+--+--"        
                Q = interpolate(Constant(0.),W_dg)
                resistance = 0.
                power = 0.
                iu = 0.
        # if (iu == nu) and (thp.em_method=="mms-nonlinear-full"):
        #     print "--+--+-- updating SAR                       --+--+--"        
        #     Q, resistance, power, Vfield = RFA_SAR(problemname, mesh, interior, boundaries, emp, T, thp)
        #     iu = 0

        # check for power < 0
        if (Q.vector().array().min()<0):
            error('source term Q < 0')


        # assemble each iteration to account for previous time step
        dte.dt = dt
        cur_time.t = t

        if thp.em_method=='iterate' or thp.em_method=='iterateRFA' or thp.em_method=="mms-nonlinear-full":
            f.assign(Q)
        D_prev_const.assign(project_axisym(D_prev,W_dg))
        b = assemble(L, tensor=b)
        A = assemble(a, tensor=A)

        for bc in bcs:
            bc.apply(A, b)
#        solve(A, T.vector(), b,
#              solver_parameters={"linear_solver": "mumps","preconditioner": "hypre_euclid"})

        solver.solve(A, T.vector(), b)

        # adapt T to ensure about 5 time steps per degree change
        nodal_T = T.vector().array()
        nodal_T_prev = T_prev.vector().array()
        T_error = N.abs(nodal_T-nodal_T_prev).max()
        print "max T err: ", T_error, "   t: ", t+dt, "   dt: ", dt, "   pow: ", power, "   imp: ", resistance

        #plot(T)
        #interactive()

        if T_error < abs(thp.Tu-thp.Tl)*0.03 and dt < dt_max and step_inc:
            # t += dt*.1
            # dt = dt*1.1
            t += dt*.1
            dt = dt*1.1
            # ensure dt_max not exceeded
            if dt > dt_max:
                # remove time step
                # dt = dt/1.1
                # t -= .1*dt
                dt = dt/1.1
                t -= .1*dt
                # replace with dt_max
                t -= dt
                t += dt_max
                dt = dt_max
            print "**************************** INCREASE STEP *********************************"
        elif T_error > abs(thp.Tu-thp.Tl)*0.3 and dt > dt_min:
#            t = t - dt*0.5
#            dt = dt*0.5
            t = t - dt*0.1
            dt = dt*0.9
            step_inc = False # stop increase immediately after decrease
            print "**************************** DECREASE STEP *********************************"
        else:
            # check that temp does not change too much
            if T_error > abs(thp.Tu-thp.Tl)*0.3:
                #print "step size too large"
                if thp.stop_on_me:
                    plot(T-T_prev)
                    interactive()
                    error("step size too large")
                else:
                    warning("step size too large")
                    print "********* TIME STEP MAX ERROR ", T_error, " ****************"
                
            # Test if output required and compute directly
            #eps = dt*0.001
#            if any(N.logical_and(t_out >= t-dt,t_out < t+eps)):
            if any(N.logical_and(t_out > t,t_out < t+dt+eps)):
                print "***************************** OUTPUT FILE **********************************"
#                dt_0=N.abs(t_out[N.logical_and(t_out>=t-dt,t_out<t+eps)]-t+dt)
                dt_0=N.abs(t_out[N.logical_and(t_out > t,t_out < t+dt+eps)]-t)
                print dt_0
                T_out.vector()[:] = T_prev.vector().array()+(T.vector().array()-T_prev.vector().array())/dt*dt_0

                # scale back to Kelvin and change name
                T_out = project_axisym(T_out+Constant(310.),W)
                T_out.rename('Temperature',T_out.label())
                file_temp << (T_out, t+dt_0[0])
                
                # rest just change name and output
                Q_out = project_axisym(f,W)
                Q_out.rename('SAR',Q_out.label())
                file_SAR << Q_out
                # V_out = project_axisym(Vfield,W)
                # V_out.rename('Voltage',V_out.label())
                # file_vfield << V_out

                D.vector()[:] = 1.-cda[1::2] # viable (also only last value not linear interp as T)
                file_cd << D

                # perf_out = project_axisym(omega,W_dg_0)
                # perf_out.rename('Perfusion',perf_out.label())
                # file_perf << perf_out


#                store_resistance[N.logical_and(t_out>=t-dt,t_out<t+eps)] = resistance
#                store_power[N.logical_and(t_out>=t-dt,t_out<t+eps)] = power
                store_resistance[N.logical_and(t_out > t,t_out < t+dt+eps)] = resistance
                store_power[N.logical_and(t_out > t,t_out < t+dt+eps)] = power
                

            # Update cell death
            print "***************************** CELL DEATH  **********************************"
            #print cda
            #print n
            #print t
            #print dt
            #print nodal_T
            #print nodal_T.min()
            #print nodal_T.max()
            if thp.cda_update:
                cda = cell_death_timestep(cda,n,t,dt,nodal_T,cell_death_parameters)
            D.vector()[:] = 1.-cda[1::2]
            D_prev.assign(D) # update cell death for perfusion
            #print cda

            # Control System
            # 
            # If Impedance greater than emp.imp_max SAR == 0
            if (resistance > emp.imp_max) and imp_on:
                imp_off_t_start = t
                imp_on = False
            #print 'imp_off_t_start', imp_off_t_start
            #print 'emp.imp_t_off', emp.imp_t_off
            if t - imp_off_t_start > emp.imp_t_off:
                imp_on = True

            #print imp_on

            # ACCEPT TIME STEP
            t += dt
            T_prev.assign(T)

            step_inc = True # allow increase after accepted step
            print "***************************** ACCEPT TIME **********************************"

    N.savetxt("%s/impedance.out" % problemname, (store_resistance, store_power), delimiter=',')   # output impedance

    return project_axisym(T+Constant(310.),W)

# cell_death_func
# 
# evaluate gradients of the ode system at a time t
#
# args:
#    y = vector of current state [alive, dead] with len(alive)==len(dead)
#    t = current time (dummy)
#    n = len(alive)
#    kb = model parameter
#    kf_ = model parameter
#    Tk = model parameter
#    T_prev_nodal = vector of corresponding temperature values len(T_pev_nodal)==len(alive)
#
# returns:
#    dydt = vector of gradients len(dydt)==len(alive)

def cell_death_func(y,t,n,kb,kf_,Tk,T_prev_nodal): #odeint
    #print y, t, n, kb,kf_,Tk,T_prev_nodal
    #dydt = N.zeros(2*n) # cell death array
    dydt = N.array(y)
    dydt[:n] = -kf_*N.exp((T_prev_nodal-273.)/Tk)*(1.-y[:n])*y[:n]+kb*(1.-y[:n]-y[n:])
    dydt[n:] = kf_*N.exp((T_prev_nodal-273.)/Tk)*(1.-y[:n])*(1.-y[:n]-y[n:])
    
    return dydt

# def cell_death_func_class(t,y,args): # ode class
#     # args = n, cdp.kb, cdp.kf_, cdp.Tk, T_prev_nodal
#     dydt = N.array(y)
#     dydt[:args[0]] = -args[2]*N.exp((args[4]-273.)/args[3])*(1.-y[:args[0]])*y[:args[0]]+args[1]*(1.-y[:args[0]]-y[args[0]:])
#     dydt[args[0]:] = args[2]*N.exp((args[4]-273.)/args[3])*(1.-y[:args[0]])*(1.-y[:args[0]]-y[args[0]:])
    
#     return dydt

def cell_death_func_class(t,y,args): # ode class
    # scaled for T' = T-310
    # args = n, cdp.kb, cdp.kf_, cdp.Tk, T_prev_nodal
    dydt = N.array(y)
    dydt[::2] = -args[2]*N.exp((args[4]+37.)/args[3])*(1.-y[::2])*y[::2]+args[1]*(1.-y[::2]-y[1::2])
    dydt[1::2] = args[2]*N.exp((args[4]+37.)/args[3])*(1.-y[::2])*(1.-y[::2]-y[1::2])
    
    return dydt

# cell_death_jac_class
#
# jacobian of the cell death model in vain hope of performance gain

def cell_death_jac_class(t,y,args):
    packed_jac = N.zeros((3,len(y)))
    # diagonal
    packed_jac[1][::2] = args[2]*N.exp((args[4]-273.)/args[3])*(2*y[::2]-1)-args[2]
    packed_jac[1][1::2] = -args[2]*N.exp((args[4]-273.)/args[3])*(1-y[::2])
    # lead (to left of diagonal)
    packed_jac[0][1::2] = -args[2]
    # trail (to right of diagonal)
    packed_jac[2][::2] = args[2]*N.exp((args[4]-273.)/args[3])*(2*y[::2]+y[1::2]-2)
    #print 'jac called'
    
    return packed_jac

# cell_death_timestep
# 
# evaluate cell death model for current time step
#
# args:
#    y0 - initial condition [alive, dead] (model state at beginning of time step)
#    n - len(alive)
#    dt - width of time step
#    T_prev_nodal - temperature to assume during time step
#    cdp - cell death parameters
#
# cell_death_timestep(N.array([0.99,0.]),1,0.,900.,338,cell_death_parameters)

def cell_death_timestep(y0,n,t,dt,T_prev_nodal,cdp):
    import scipy.integrate as spint
    # ode int
    # time = N.array([t,t+dt])
    # yt = spint.odeint(cell_death_func,y0,time,(n,cdp.kb,cdp.kf_,cdp.Tk,T_prev_nodal),atol=1e-4,rtol=1e-4,mxords=4)
    # cda = yt[1]
    # ode class
    step_int = spint.ode(cell_death_func_class)
    step_int.set_integrator("vode",method="bdf",nsteps=1e4)

    # step_int = spint.ode(cell_death_func_class,cell_death_jac_class)
    # step_int.set_integrator("vode",method="bdf",nsteps=1e4,lband=1,with_jacobian=True)
    # step_int.set_jac_params([n, cdp.kb, cdp.kf_, cdp.Tk, T_prev_nodal])

    step_int.set_f_params([n, cdp.kb, cdp.kf_, cdp.Tk, T_prev_nodal])
    step_int.set_initial_value(y0,0.)
    step_int.integrate(dt)
    print step_int.successful()
    if not step_int.successful():
        error("cell death step solve failed")
    cda = step_int.y

    return cda

# RFA_SAR
#
# compute the SAR according to the Laplace equation in axisymmetric
# cyclindrical coordinates
#
# args:
#    

def RFA_SAR(problemname, mesh, interior, boundaries, emp, T, thp):
    
    # public function pre-amble
    ensure_dir(problemname) # check directory exists for results

    print "--+--+-- compute RFA SAR --+--+--"        

    # set solver params in advance
    solver = KrylovSolver("cg", "hypre_euclid")
    #solver.parameters["absolute_tolerance"] = 1E-7
    #solver.parameters["relative_tolerance"] = 1E-4
    #solver.parameters["maximum_iterations"] = 1000
    #set_log_level(DEBUG)

    # define a restriction for V calculation
    restriction = Restriction(interior,emp.restrict_mesh) # restrict em calc
    W = FunctionSpace(restriction, 'CG', order)
    W_dg0 = FunctionSpace(restriction, 'DG', 0) # attempt to stabilise sigma
    W_dg = FunctionSpace(restriction, 'DG', order) # capture SAR shape and discontinuity

    # interpolate T onto piecewise constant to stabilise conductivity
    T_p = project_axisym(T,W_dg0)

    # set measure
    dss = ds[boundaries]

    # Get the r and z components
    polar = W.cell().x
    r = polar[0]
    z = polar[1]
    
    # allocate boundary conditions to line regions defined in boundary mesh function
    # symmetry and insulating are natural BCs
    bcs = []
    for index in emp.zero_V:
        bcs.append(DirichletBC(W, emp.V0, boundaries, index))
    
    for index in emp.active_tip:
        bcs.append(DirichletBC(W, emp.Vprobe, boundaries, index))

    # define variational problem
    V = TrialFunction(W)
    U = TestFunction(W)

    # define electrical conductivity
    sigma=emp.cond
    if emp.cond_model=='linear':
        sigma=emp.cond + emp.cond_rate*(T_p)
    elif emp.cond_model=='nonlinear':
        sigma=(conditional(le(T_p,thp.Tu), emp.cond + emp.cond_rate*(T_p), 0.) +
            conditional(And(gt(T_p,thp.Tu),le(T_p,thp.Tu+5.)), (emp.cond_vap - (emp.cond + emp.cond_rate*(thp.Tu)))/5.*(T_p-thp.Tu) + (emp.cond + emp.cond_rate*(thp.Tu)), 0.) +
            conditional(gt(T_p,thp.Tu+5.), emp.cond_vap, 0.))

#    File("%s/sigma.pvd" % problemname) << project_axisym(sigma,W_dg0)

    a = sigma*inner(nabla_grad(V), nabla_grad(U))*r*dx
    L = emp.mms_source*U*r*dx

    V = Function(W)

    solve(a == L, V, bcs)

    U = TestFunction(W_dg)
    SAR = TrialFunction(W_dg)
    a = SAR*U*r*dx
    L = U*sigma*inner(nabla_grad(V), nabla_grad(V))*r*dx

    SAR = Function(W_dg)
   
    solve(a == L, SAR)

    # SAR = project_axisym(sigma*inner(nabla_grad(V), nabla_grad(V))*r,W_dg)

    # plot(V)
    # plot(SAR)
    # interactive()

    # compute impedance for control system as purely electrical
    # need to use assemble
    # integrate current density over probe active tip
    # line element should be ds without adaption

    power = assemble(SAR*r*dx)*2*N.pi
    resistance = (emp.Vprobe)**2/power
#    print "power: ", power
#    print "resistance according to Kroger & elmer: ", (emp.Vprobe)**2/power

    return SAR, resistance, power, V

# define a projection function for axisymetric to avoid mistakes

def project_axisym(func,space):
    polar = space.cell().x
    r = polar[0]
    z = polar[1]

    w = TrialFunction(space)
    v = TestFunction(space)

    a = inner(w,v)*r*dx
    L = inner(func, v)*r*dx

    pfunc = Function(space)
    solve(a == L, pfunc)

    return pfunc
