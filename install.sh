#/bin/bash

for f in tests/*.py; do echo "processing $f"; filename=$(basename $f ".py"); gmsh -2 "mesh/$filename.geo"; dolfin-convert "mesh/$filename.msh" "mesh/$filename.xml"; done
