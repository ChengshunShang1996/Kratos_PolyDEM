#  Kratos Multi-Physics - DEM Application
#
#  License:       BSD License
#                 Kratos default license: kratos/license.txt
#
#  Main authors:  Chengshun Shang (cshang@cimne.upc.edu)
#                 
#

import numpy as np
import weakref
import shutil
import vtk
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

class PolyhedronVtkOutput():
    def __init__(self, main_path, problem_name, spheres_model_part, contact_model_part, rigid_face_model_part, polyhedron_model_part, DEM_parameters):
        self.problem_name = problem_name
        self.DEM_parameters = DEM_parameters

        # Reading Post options from DEM_parameters
        self.PostDisplacement = self.DEM_parameters["PostDisplacement"].GetBool()
        self.PostVelocity = self.DEM_parameters["PostVelocity"].GetBool()  

        self.spheres_model_part = weakref.proxy(spheres_model_part)
        self.contact_model_part = weakref.proxy(contact_model_part)
        self.rigid_face_model_part = weakref.proxy(rigid_face_model_part)
        self.polyhedron_model_part = weakref.proxy(polyhedron_model_part)

        self.vtk_post_path_directory = os.path.join(main_path, problem_name + "_Post_POLY_VTK_Files")
        shutil.rmtree(self.vtk_post_path_directory, ignore_errors=True)
        os.makedirs(str(self.vtk_post_path_directory))

        self.counter = 0

    def DataPreparation(self):

        number_of_nodes = self.polyhedron_model_part.NumberOfNodes(0)

        self.polygon_centers = np.zeros((number_of_nodes, 3))
        self.polygon_origins = []

        i = 0
        for node in self.polyhedron_model_part.Nodes:
            self.polygon_centers[i] = [node.X, node.Y, node.Z]
            i += 1
        
        self.polygon_origins.append(np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]))
        self.polygon_origins.append(np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]))
    
    def WriteResults(self, time):

        self.DataPreparation()
        
        points_vtk = vtk.vtkPoints()
        cells_vtk = vtk.vtkCellArray()

        polyhedron_points = []
        i = 0
        for polygon_center in self.polygon_centers:
            polygon = polygon_center + self.polygon_origins[i]
            for poly in polygon:
                polyhedron_points.append(poly)
            i += 1

        polyhedron_faces = []
        faces_origin = "(0,1,2,3) (1,5,6,2) (5,6,7,4) (4,7,3,0) (3,2,6,7) (0,1,5,4)"
        faces = faces_origin.strip().split()
        
        offset = 0
        for face in faces:
            vertices = list(map(int, face.strip('()').split(',')))
            adjusted_vertices = [v + offset for v in vertices]
            polyhedron_faces.append(adjusted_vertices)
        
        offset += 8

        for face in faces:
            vertices = list(map(int, face.strip('()').split(',')))
            adjusted_vertices = [v + offset for v in vertices]
            polyhedron_faces.append(adjusted_vertices)
        
        for point_coord in polyhedron_points:
            points_vtk.InsertNextPoint(point_coord)

        for face in polyhedron_faces:
            cell = vtk.vtkPolygon()
            cell.GetPointIds().SetNumberOfIds(len(face))
            for j, point_id in enumerate(face):
                cell.GetPointIds().SetId(j, point_id)
            cells_vtk.InsertNextCell(cell)

        grid = vtk.vtkUnstructuredGrid()
        grid.SetPoints(points_vtk)
        grid.SetCells(vtk.VTK_POLYGON, cells_vtk)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        file_path = os.path.join(self.vtk_post_path_directory, "polyhedron_{}.vtu".format(time))
        writer.SetFileName(file_path)
        writer.SetInputData(grid)
        writer.Write()
    