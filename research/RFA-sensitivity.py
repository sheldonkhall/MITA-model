# RFA-sensitivity
#
# input file for running the sensitivity analysis
#
# This file repeatedly calls the MITA model using a table of input parameters
# that explore the parameter space. The analysis of the results of these
# experiments is what will produce the sensitivity analysis ranking
# parameters in order of their importance.
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
from scipy import interpolate
import time as tm
import os
import scipy.io as sio

problemname = "RFA-sensitivity"

# define general parameters
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

# Load geometry
mesh = Mesh("mesh/%s.xml" % problemname)
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % problemname)
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % problemname)

# set thermal parameters
thermal_parameters.rho_c_t = 4.0e6 # normal tissue
thermal_parameters.rho_c_v = 0.6e6 # values for vapourised tissue
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

# solver options
dt_min = .00001 # absolute step size minimum (0.0001 good)
dt_max = .01 # absolute step size maximum
# tmax = 600 # maximum time (s)
# t_out = np.linspace(1,tmax,600) # numpy vector of times at which to save to disk
tmax = 25 # maximum time (s)
t_out = np.linspace(0,tmax,1000) # numpy vector of times at which to save to disk
dt = .00001 # time step (s)
thermal_parameters.em_method = 'iterateRFA'
thermal_parameters.k_model = 'linear_limited'
thermal_parameters.perf_model = 'stop'
print "need to check version with nonlinear perfusion, may not be updating"

# CHANGE
thermal_parameters.stop_on_me = False


set_log_active(False) # switch off fenics messages

# pick experiment to run
experiment = 'experiment-1'

# define function for parameter scaling
def scale_param(upper,lower,value):
    change = upper - lower
    change = change*value
    change = lower + change
    return change

# define parameter ranges
rho_c_t = 3.7e6
rho_c_t_ = 4.3e6
rho_c_v = .44e6
rho_c_v_ = .8e6
k = .46
k_ = .57
dk = 0.
dk_ = .0033
cond = .14
cond_ = .28
cond_rate = 1
cond_rate_ = 2
omega = 34020
omega_ = 68040
Cliq = .71
Cliq_ = .76
dT = 1.
dT_ = 10.
pc = np.linspace(0,1,5)
#kf_ = .8e-3
#kf__ = 9.1e-3
kf = interpolate.interp1d(pc,np.array([.0008, .00262, .00352, .00454, .00907]))
#kb = .25e-3
#kb_ = 19.e-3
kb = interpolate.interp1d(pc,np.array([.00025, .00574, .00846, .0108, .0192]))
#Tk = 25.
#Tk_ = 64.
Tk = interpolate.interp1d(pc,np.array([24.6, 36.7, 41.6, 46.3, 63.5]))
cond_vap = 0.01
cond_vap_ = 0.0001
p_stop = 0.7
p_stop_ = 0.9



# load experiment
expdes = sio.loadmat(problemname + '/' + experiment + '.mat')
parameters = expdes['A']

# change to specific directory for experiment
file_prefix = tm.strftime('%Y%m%d-%H%M%S')
if not os.path.exists('%s' % problemname):
        os.makedirs('%s' % problemname)
problemname = problemname + '/' + file_prefix + '-' + experiment
os.makedirs('%s' % (problemname))
logfile = open('%s/log.txt' % problemname, 'w')
logfile.write('Beginning ' + experiment + '\n')

# loop through runs in experiment
#for run in range(parameters.shape[0]):
for run in [0]:
    
    # make dir for run
    problemnametemp = problemname + '/run' + str(run+1)
    os.makedirs('%s' % (problemnametemp))
    
    # start time of run
    logfile.write('Beginning run ' + str(run+1) + ' at ' + tm.strftime('%H:%M:%S') + '\n')

    # set parameters
    thermal_parameters.rho_c_t = scale_param(rho_c_t_, rho_c_t, parameters[run,0])
    logfile.write('rho_c_t: ' + str(scale_param(rho_c_t_, rho_c_t, parameters[run,0])) + '\n')
    thermal_parameters.rho_c_v = scale_param(rho_c_v_, rho_c_v, parameters[run,1])
    logfile.write('rho_c_v: ' + str(scale_param(rho_c_v_, rho_c_v, parameters[run,1])) + '\n')
    thermal_parameters.k = scale_param(k_, k, parameters[run,2])
    logfile.write('k: ' + str(scale_param(k_, k, parameters[run,2])) + '\n')
    thermal_parameters.dk = scale_param(dk_, dk, parameters[run,3])
    logfile.write('dk: ' + str(scale_param(dk_, dk, parameters[run,3])) + '\n')
    EM_parameters.cond = scale_param(cond_, cond, parameters[run,4])
    logfile.write('cond: ' + str(scale_param(cond_, cond, parameters[run,4])) + '\n')
    EM_parameters.cond_rate = scale_param(cond_rate_, cond_rate, parameters[run,5])
    logfile.write('cond_rate: ' + str(scale_param(cond_rate_, cond_rate, parameters[run,5])) + '\n')
    thermal_parameters.omega = scale_param(omega_, omega, parameters[run,6])
    logfile.write('omega: ' + str(scale_param(omega_, omega, parameters[run,6])) + '\n')
    thermal_parameters.Cliq = scale_param(Cliq_, Cliq, parameters[run,7])
    logfile.write('Cliq: ' + str(scale_param(Cliq_, Cliq, parameters[run,7])) + '\n')
    ###  problem with Tu and Tl being overwritten in enthalpy function
#    thermal_parameters.Tl = thermal_parameters.Tu - scale_param(dT_, dT, parameters[run,8])
    thermal_parameters.Tu = 373
    thermal_parameters.Tl = (373 - scale_param(dT_, dT, parameters[run,8]))
    logfile.write('dT: ' + str(scale_param(dT_, dT, parameters[run,8])) + '\n')
    # cell_death_parameters.kf_ = scale_param(kf__, kf_, parameters[run,9])
    # logfile.write('kf_: ' + str(scale_param(kf__, kf_, parameters[run,9])) + '\n')
    # cell_death_parameters.kb = scale_param(kb_, kb, parameters[run,10])
    # logfile.write('kb: ' + str(scale_param(kb_, kb, parameters[run,10])) + '\n')
    # cell_death_parameters.Tk = scale_param(Tk_, Tk, parameters[run,11])
    # logfile.write('Tk: ' + str(scale_param(Tk_, Tk, parameters[run,11])) + '\n')
    cell_death_parameters.kf_ = kf(parameters[run,9])
    logfile.write('kf_: ' + str(cell_death_parameters.kf_) + '\n')
    cell_death_parameters.kb = kb(parameters[run,9])
    logfile.write('kb: ' + str(cell_death_parameters.kb) + '\n')
    cell_death_parameters.Tk = Tk(parameters[run,9])
    logfile.write('Tk: ' + str(cell_death_parameters.Tk) + '\n')
    EM_parameters.cond_vap = scale_param(cond_vap_, cond_vap, parameters[run,10])*EM_parameters.cond
    logfile.write('cond_vap: ' + str(scale_param(cond_vap_, cond_vap, parameters[run,10])*EM_parameters.cond) + '\n')
    thermal_parameters.p_stop = scale_param(p_stop_, p_stop, parameters[run,11])
    logfile.write('p_stop: ' + str(scale_param(p_stop_, p_stop, parameters[run,11])) + '\n')
    
    # force write output to logfile
    logfile.flush()

    # call model
    T = compute_enthalpy_nl(mesh, interior, boundaries, problemnametemp, dt, tmax, dt_min, dt_max, t_out, thermal_parameters, EM_parameters)

    # end time
    logfile.write('Ended at: ' + tm.strftime('%H:%M:%S') + '\n')

logfile.write('--- End of Experiment ---')
