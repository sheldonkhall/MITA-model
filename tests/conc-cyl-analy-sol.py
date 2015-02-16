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

#import numpy as N
from scipy import constants as S
from scipy import special
import numpy as N
from dolfin import *
import sys

#
# - Some constants
#

pi = S.pi
e = S.e
c0 = S.c
mu_0 = S.mu_0
eps_0 = S.epsilon_0

#
# - Define Problem here
#

problemname = "conc-cyl-test"
set_log_level(ERROR) # remove warnings for tests

# Define uniform problemname parameters
freq = 7.0e7 # Probe frequency
om = 2 * pi * freq # Probe angular frequency

# Define domain-specific problemname parameters
# (eps_r, sig, mu_r)
tissue = [[0.0 for x in range(2)] for y in range(12)]
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

# ([dummy], inner annulus (n), outer annulus (m))
n = 11
m = 2
eps_r_by_subdomain = [-1, tissue[n][0], tissue[m][0]]
sigma_by_subdomain = [-1, tissue[n][1], tissue[m][1]]
mu_r_by_subdomain = [-1, tissue[n][2], tissue[m][2]]



# compute analytic solution

# load variables

ea1 = N.array(0.0,dtype='complex')
ea1 = eps_r_by_subdomain[1]*eps_0 + 1j*sigma_by_subdomain[1]/om
ea2 = N.array(0.0,dtype='complex')
ea2 = eps_r_by_subdomain[2]*eps_0 + 1j*sigma_by_subdomain[2]/om
k1 = om*N.sqrt(mu_0*ea1)
k2 = om*N.sqrt(mu_0*ea2)
a = 0.12
R = 0.25

# construct system of equations
## A = N.matrix(N.zeros([3,3],dtype='complex'))
## b = N.matrix(N.zeros([3,1],dtype='complex'))
A = N.array(N.zeros([3,3],dtype='complex'))
b = N.array(N.zeros([3,1],dtype='complex'))

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

print A

x = N.linalg.solve(A,b)

print x

print N.dot(A,x)
