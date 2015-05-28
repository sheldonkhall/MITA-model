#!/usr/bin/env ipython
#
# open vtu file from fenics/elmer and compute the area of lesion in
# axisymmetric calculations.
#
# by Sheldon Hall
 
from vtk import *
import sys
import scipy.io

def IsoAreaCompute( file_name, thresh, scField ):
    # char file_name - vtu file containing solution field
    # real thresh - threshold for lesion
    # char scField - name of active scalar field

    #scField = 'dead'
    #scField = 'f_64'

    # Create renderer and add actors to it
    renderer = vtkRenderer()
    renderer.SetBackground( 1.0, 1.0 , 1.0 )

    # Create render window
    window = vtkRenderWindow()
    window.AddRenderer( renderer )
    window.SetSize( 500, 500 )

    # Create interactor
    interactor = vtkRenderWindowInteractor()
    interactor.SetRenderWindow( window )

    # Read the source file.
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update() # Needed because of GetScalarRange
    output = reader.GetOutput()

    # convert unstructured to poly
    conv = vtkGeometryFilter()
    conv.SetInput( output )
    conv.Update()
    poly_data = vtkPolyData()
    poly_data = conv.GetOutput()

    # select array of interest
#    print output.GetPointData().GetScalars().GetName()
#    poly_data.GetPointData().SetActiveScalars( scField )
    poly_data.GetPointData().SetActiveScalars( output.GetPointData().GetScalars().GetName() )

    # clip region using threshold value
    cut = vtkClipPolyData()
    cut.SetInput( poly_data )
    cut.SetValue( thresh )
    cut.InsideOutOn() # compute area of scalar field under thresh
    # cut.InsideOutOff() # compute area of scalar field over thresh
    
    cut.Update()
    mesh = cut.GetOutput()

    # compute area from polydata
    area = vtkMassProperties()
    area.SetInput( mesh )
    area.Update()
    lesion = area.GetSurfaceArea()

    scipy.io.savemat('lesion_area.mat', mdict={'lesion_area': lesion})

    # extract temperature range
    #scalar_range = mesh.GetScalarRange()

    # Create the mapper for isovolume
    #meshMapper = vtkDataSetMapper()
    #meshMapper.SetInput( mesh )
    #meshMapper.SetScalarRange( scalar_range )
    #meshMapper.SetScalarRange( (0., scalar_range[1]) )

    # Create wireframe actor for submesh
    #meshActor = vtkActor()
    #meshActor.SetMapper( meshMapper )
    #meshActor.GetProperty().SetRepresentationToWireframe()
    #renderer.AddActor( meshActor )

    # Start interaction
    #window.Render()
    #interactor.Start()
    #del interactor
    #del window

    return lesion

if __name__ == '__main__':
    a = sys.argv[1]
    c = sys.argv[2]
    d = sys.argv[3]
    b = IsoAreaCompute(a,float(c),d)
    sys.stdout.write(str(b))
    sys.stdout.write('\n')
    sys.stdout.flush()
    sys.exit(0)
