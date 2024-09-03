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
        self.polyhedron_faces = []
        self.polygon_velocity = np.zeros((number_of_nodes, 3))
        self.polygon_displacement = np.zeros((number_of_nodes, 3))

        i = 0
        for node in self.polyhedron_model_part.Nodes:
            self.polygon_centers[i] = [node.X, node.Y, node.Z]
            i += 1

        for element in self.polyhedron_model_part.Elements:
            vertices = element.GetListOfVertices()
            temp_vertices_list = []
            for vertex in vertices:
                temp_vertices_list.append([vertex[0], vertex[1], vertex[2]])
            self.polygon_origins.append(temp_vertices_list)

        for element in self.polyhedron_model_part.Elements:
            faces = element.GetListOfFaces()
            temp_faces_list = []
            for face in faces:
                temp_face_sub = []
                for face_sub in face:
                    temp_face_sub.append(face_sub)
                temp_faces_list.append(temp_face_sub)
            self.polyhedron_faces.append(temp_faces_list)

        if self.PostVelocity:
            i = 0
            for node in self.polyhedron_model_part.Nodes:
                self.polygon_velocity[i] = [node.GetSolutionStepValue(VELOCITY_X), node.GetSolutionStepValue(VELOCITY_Y), node.GetSolutionStepValue(VELOCITY_Z)]
                i += 1

        if self.PostDisplacement:
            i = 0
            for node in self.polyhedron_model_part.Nodes:
                self.polygon_displacement[i] = [node.GetSolutionStepValue(DISPLACEMENT_X), node.GetSolutionStepValue(DISPLACEMENT_Y), node.GetSolutionStepValue(DISPLACEMENT_Z)]
                i += 1
    
    def WriteResults(self, poly_output_cnt):

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
        i = 0
        offset = 0
        for polyhedron_face in self.polyhedron_faces:
            if i > 0:
                offset += len(self.polygon_origins[i-1])
            for face in polyhedron_face:
                adjusted_face = [v + offset for v in face]
                polyhedron_faces.append(adjusted_face)
            i += 1
        
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

        if self.PostVelocity:
            velocity_array = vtk.vtkDoubleArray()
            velocity_array.SetName("Velocity")
            velocity_array.SetNumberOfComponents(3)

            i = 0
            for ii in range(len(self.polygon_centers)):
                for jj in range(len(self.polygon_origins[i])):
                    velocity_array.InsertNextTuple(self.polygon_velocity[i])
                i += 1

            grid.GetPointData().AddArray(velocity_array)

        if self.PostDisplacement:
            displacement_array = vtk.vtkDoubleArray()
            displacement_array.SetName("Displacement")
            displacement_array.SetNumberOfComponents(3)

            i = 0
            for ii in range(len(self.polygon_centers)):
                for jj in range(len(self.polygon_origins[i])):
                    displacement_array.InsertNextTuple(self.polygon_displacement[i])
                i += 1

            grid.GetPointData().AddArray(displacement_array)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        file_path = os.path.join(self.vtk_post_path_directory, "polyhedron_{}.vtu".format(poly_output_cnt))
        writer.SetFileName(file_path)
        writer.SetInputData(grid)
        writer.Write()
    