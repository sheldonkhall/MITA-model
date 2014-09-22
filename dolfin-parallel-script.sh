#!/usr/bin/ipython

# take a fenics meshfunction from gmsh in xml format and write over
# in new xml format - workaround for bug in parallel applications

import sys
from dolfin import *

mesh = Mesh("mesh/%s.xml" % sys.argv[1])
boundaries = MeshFunction("size_t", mesh, "mesh/%s_facet_region.xml" % sys.argv[1])
interior = MeshFunction("size_t", mesh, "mesh/%s_physical_region.xml" % sys.argv[1])

File("mesh/%s_facet_region.xml" % sys.argv[1]) << boundaries
File("mesh/%s_physical_region.xml" % sys.argv[1]) << interior
