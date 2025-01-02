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
        self.PostAngularVelocity = self.DEM_parameters["PostAngularVelocity"].GetBool()  
        self.PostTotalForces = self.DEM_parameters["PostTotalForces"].GetBool() 

        self.spheres_model_part = weakref.proxy(spheres_model_part)
        self.contact_model_part = weakref.proxy(contact_model_part)
        self.rigid_face_model_part = weakref.proxy(rigid_face_model_part)
        self.polyhedron_model_part = weakref.proxy(polyhedron_model_part)

        self.vtk_post_path_directory = os.path.join(main_path, problem_name + "_Post_POLY_VTK_Files")
        shutil.rmtree(self.vtk_post_path_directory, ignore_errors=True)
        os.makedirs(str(self.vtk_post_path_directory))

        self.counter = 0

    def DataPreparation(self):

        self.find_smp_p = False
        for smp in self.polyhedron_model_part.SubModelParts:
            if smp.Name == "DEMParts_Body":
                self.polyhedron_model_part_p = smp
                self.find_smp_p = True

        self.find_smp_w = False
        for smp in self.polyhedron_model_part.SubModelParts:
            if smp.Name == "DEMWall":
                self.polyhedron_model_part_w = smp
                self.find_smp_w = True
        
        if self.find_smp_p and self.find_smp_w:
            self.DataPreparationParticles()
            self.DataPreparationWalls()
        else:
            self.DataPreparationAll()

    def DataPreparationAll(self):

        number_of_nodes = self.polyhedron_model_part.NumberOfNodes(0)

        self.polygon_centers = np.zeros((number_of_nodes, 3))
        self.polygon_origins = []
        self.polyhedron_faces = []
        self.polygon_velocity = np.zeros((number_of_nodes, 3))
        self.polygon_angular_velocity = np.zeros((number_of_nodes, 3))
        self.polygon_displacement = np.zeros((number_of_nodes, 3))
        self.polygon_totalforces = np.zeros((number_of_nodes, 3))
        self.polygon_id = np.zeros((number_of_nodes, 1))

        i = 0
        for node in self.polyhedron_model_part.Nodes:
            self.polygon_centers[i] = [node.X, node.Y, node.Z]
            i += 1
        
        i = 0
        for node in self.polyhedron_model_part.Nodes:
            self.polygon_id[i] = node.Id
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
        
        if self.PostAngularVelocity:
            i = 0
            for node in self.polyhedron_model_part.Nodes:
                self.polygon_angular_velocity[i] = [node.GetSolutionStepValue(ANGULAR_VELOCITY_X), node.GetSolutionStepValue(ANGULAR_VELOCITY_Y), node.GetSolutionStepValue(ANGULAR_VELOCITY_Z)]
                i += 1

        if self.PostDisplacement:
            i = 0
            for node in self.polyhedron_model_part.Nodes:
                self.polygon_displacement[i] = [node.GetSolutionStepValue(DISPLACEMENT_X), node.GetSolutionStepValue(DISPLACEMENT_Y), node.GetSolutionStepValue(DISPLACEMENT_Z)]
                i += 1

        if self.PostTotalForces:
            i = 0
            for node in self.polyhedron_model_part.Nodes:
                self.polygon_totalforces[i] = [node.GetSolutionStepValue(TOTAL_FORCES)[0], node.GetSolutionStepValue(TOTAL_FORCES)[1], node.GetSolutionStepValue(TOTAL_FORCES)[2]]
                i += 1

    def DataPreparationParticles(self):
        
        number_of_nodes_p = self.polyhedron_model_part_p.NumberOfNodes(0)

        self.polygon_centers_p = np.zeros((number_of_nodes_p, 3))
        self.polygon_origins_p = []
        self.polyhedron_faces_p = []
        self.polygon_velocity_p = np.zeros((number_of_nodes_p, 3))
        self.polygon_angular_velocity_p = np.zeros((number_of_nodes_p, 3))
        self.polygon_displacement_p = np.zeros((number_of_nodes_p, 3))
        self.polygon_totalforces_p = np.zeros((number_of_nodes_p, 3))
        self.polygon_id_p = np.zeros((number_of_nodes_p, 1))

        i = 0
        for node in self.polyhedron_model_part_p.Nodes:
            self.polygon_centers_p[i] = [node.X, node.Y, node.Z]
            i += 1

        i = 0
        for node in self.polyhedron_model_part_p.Nodes:
            self.polygon_id_p[i] = node.Id
            i += 1

        for element in self.polyhedron_model_part_p.Elements:
            vertices = element.GetListOfVertices()
            temp_vertices_list = []
            for vertex in vertices:
                temp_vertices_list.append([vertex[0], vertex[1], vertex[2]])
            self.polygon_origins_p.append(temp_vertices_list)

        for element in self.polyhedron_model_part_p.Elements:
            faces = element.GetListOfFaces()
            temp_faces_list = []
            for face in faces:
                temp_face_sub = []
                for face_sub in face:
                    temp_face_sub.append(face_sub)
                temp_faces_list.append(temp_face_sub)
            self.polyhedron_faces_p.append(temp_faces_list)

        if self.PostVelocity:
            i = 0
            for node in self.polyhedron_model_part_p.Nodes:
                self.polygon_velocity_p[i] = [node.GetSolutionStepValue(VELOCITY_X), node.GetSolutionStepValue(VELOCITY_Y), node.GetSolutionStepValue(VELOCITY_Z)]
                i += 1

        if self.PostAngularVelocity:
            i = 0
            for node in self.polyhedron_model_part_p.Nodes:
                self.polygon_angular_velocity_p[i] = [node.GetSolutionStepValue(ANGULAR_VELOCITY_X), node.GetSolutionStepValue(ANGULAR_VELOCITY_Y), node.GetSolutionStepValue(ANGULAR_VELOCITY_Z)]
                i += 1

        if self.PostDisplacement:
            i = 0
            for node in self.polyhedron_model_part_p.Nodes:
                self.polygon_displacement_p[i] = [node.GetSolutionStepValue(DISPLACEMENT_X), node.GetSolutionStepValue(DISPLACEMENT_Y), node.GetSolutionStepValue(DISPLACEMENT_Z)]
                i += 1

        if self.PostTotalForces:
            i = 0
            for node in self.polyhedron_model_part_p.Nodes:
                self.polygon_totalforces_p[i] = [node.GetSolutionStepValue(TOTAL_FORCES)[0], node.GetSolutionStepValue(TOTAL_FORCES)[1], node.GetSolutionStepValue(TOTAL_FORCES)[2]]
                i += 1

    def DataPreparationWalls(self):
        
        number_of_nodes_w = self.polyhedron_model_part_w.NumberOfNodes(0)

        self.polygon_centers_w = np.zeros((number_of_nodes_w, 3))
        self.polygon_origins_w = []
        self.polyhedron_faces_w = []
        self.polygon_velocity_w = np.zeros((number_of_nodes_w, 3))
        self.polygon_angular_velocity_w = np.zeros((number_of_nodes_w, 3))
        self.polygon_displacement_w = np.zeros((number_of_nodes_w, 3))
        self.polygon_totalforces_w = np.zeros((number_of_nodes_w, 3))
        self.polygon_id_w = np.zeros((number_of_nodes_w, 1))

        i = 0
        for node in self.polyhedron_model_part_w.Nodes:
            self.polygon_centers_w[i] = [node.X, node.Y, node.Z]
            i += 1

        i = 0
        for node in self.polyhedron_model_part_w.Nodes:
            self.polygon_id_w[i] = node.Id
            i += 1

        for element in self.polyhedron_model_part_w.Elements:
            vertices = element.GetListOfVertices()
            temp_vertices_list = []
            for vertex in vertices:
                temp_vertices_list.append([vertex[0], vertex[1], vertex[2]])
            self.polygon_origins_w.append(temp_vertices_list)
        
        for element in self.polyhedron_model_part_w.Elements:
            faces = element.GetListOfFaces()
            temp_faces_list = []
            for face in faces:
                temp_face_sub = []
                for face_sub in face:
                    temp_face_sub.append(face_sub)
                temp_faces_list.append(temp_face_sub)
            self.polyhedron_faces_w.append(temp_faces_list)
        
        if self.PostVelocity:
            i = 0
            for node in self.polyhedron_model_part_w.Nodes:
                self.polygon_velocity_w[i] = [node.GetSolutionStepValue(VELOCITY_X), node.GetSolutionStepValue(VELOCITY_Y), node.GetSolutionStepValue(VELOCITY_Z)]
                i += 1

        if self.PostAngularVelocity:
            i = 0
            for node in self.polyhedron_model_part_w.Nodes:
                self.polygon_angular_velocity_w[i] = [node.GetSolutionStepValue(ANGULAR_VELOCITY_X), node.GetSolutionStepValue(ANGULAR_VELOCITY_Y), node.GetSolutionStepValue(ANGULAR_VELOCITY_Z)]
                i += 1

        if self.PostDisplacement:
            i = 0
            for node in self.polyhedron_model_part_w.Nodes:
                self.polygon_displacement_w[i] = [node.GetSolutionStepValue(DISPLACEMENT_X), node.GetSolutionStepValue(DISPLACEMENT_Y), node.GetSolutionStepValue(DISPLACEMENT_Z)]
                i += 1 

        if self.PostTotalForces:
            i = 0
            for node in self.polyhedron_model_part_w.Nodes:
                self.polygon_totalforces_w[i] = [node.GetSolutionStepValue(TOTAL_FORCES)[0], node.GetSolutionStepValue(TOTAL_FORCES)[1], node.GetSolutionStepValue(TOTAL_FORCES)[2]]
                i += 1    
    
    def WriteOutAll(self, poly_output_cnt):

        polyhedron_points = []
        i = 0
        for polygon_center in self.polygon_centers:
            polygon = polygon_center + self.polygon_origins[i]
            for poly in polygon:
                polyhedron_points.append(poly)
            i += 1

        polyhedron_faces_list = []
        polyhedron_faces = []
        i = 0
        offset = 0
        for polyhedron_face in self.polyhedron_faces:
            if i > 0:
                offset += len(self.polygon_origins[i-1])
            for face in polyhedron_face:
                adjusted_face = [v + offset for v in face]
                polyhedron_faces.append(adjusted_face)
            polyhedron_faces_list.append(polyhedron_faces)
            polyhedron_faces = []
            i += 1
        
        points_vtk = vtk.vtkPoints()
        for point_coord in polyhedron_points:
            points_vtk.InsertNextPoint(point_coord)
        
        grid = vtk.vtkUnstructuredGrid()

        for poly_id, polyhedron_faces in enumerate(polyhedron_faces_list):

            faces_id_list = vtk.vtkIdList()
            faces_id_list.InsertNextId(len(polyhedron_faces))  # Number of faces
            for face in polyhedron_faces:
                faces_id_list.InsertNextId(len(face))  # Number of vertices in this face
                for idx in face:
                    faces_id_list.InsertNextId(idx)

            grid.InsertNextCell(vtk.VTK_POLYHEDRON, faces_id_list)

        grid.SetPoints(points_vtk)

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

        if self.PostAngularVelocity:
            angular_velocity_array = vtk.vtkDoubleArray()
            angular_velocity_array.SetName("Angular Velocity")
            angular_velocity_array.SetNumberOfComponents(3)

            i = 0
            for ii in range(len(self.polygon_centers)):
                for jj in range(len(self.polygon_origins[i])):
                    angular_velocity_array.InsertNextTuple(self.polygon_angular_velocity[i])
                i += 1

            grid.GetPointData().AddArray(angular_velocity_array)

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

        if self.PostTotalForces:
            totalforces_array = vtk.vtkDoubleArray()
            totalforces_array.SetName("Total Forces")
            totalforces_array.SetNumberOfComponents(3)

            i = 0
            for ii in range(len(self.polygon_centers)):
                for jj in range(len(self.polygon_origins[i])):
                    totalforces_array.InsertNextTuple(self.polygon_totalforces[i])
                i += 1

            grid.GetPointData().AddArray(totalforces_array)

        polygon_id_array = vtk.vtkIntArray()
        polygon_id_array.SetName("Polygon Id")
        polygon_id_array.SetNumberOfComponents(1)

        i = 0
        for ii in range(len(self.polygon_centers)):
            for jj in range(len(self.polygon_origins[i])):
                polygon_id_array.InsertNextTuple(self.polygon_id[i])
            i += 1

        grid.GetPointData().AddArray(polygon_id_array)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        file_path = os.path.join(self.vtk_post_path_directory, "polyhedron_{}.vtu".format(poly_output_cnt))
        writer.SetFileName(file_path)
        writer.SetInputData(grid)
        writer.Write()

    def WriteOutParticles(self, poly_output_cnt):

        polyhedron_points = []
        i = 0
        for polygon_center in self.polygon_centers_p:
            polygon = polygon_center + self.polygon_origins_p[i]
            for poly in polygon:
                polyhedron_points.append(poly)
            i += 1

        polyhedron_faces_list = []
        polyhedron_faces = []
        i = 0
        offset = 0
        for polyhedron_face in self.polyhedron_faces_p:
            if i > 0:
                offset += len(self.polygon_origins_p[i-1])
            for face in polyhedron_face:
                adjusted_face = [v + offset for v in face]
                polyhedron_faces.append(adjusted_face)
            polyhedron_faces_list.append(polyhedron_faces)
            polyhedron_faces = []
            i += 1
        
        points_vtk = vtk.vtkPoints()
        for point_coord in polyhedron_points:
            points_vtk.InsertNextPoint(point_coord)
        
        grid = vtk.vtkUnstructuredGrid()

        for poly_id, polyhedron_faces in enumerate(polyhedron_faces_list):

            faces_id_list = vtk.vtkIdList()
            faces_id_list.InsertNextId(len(polyhedron_faces))  # Number of faces
            for face in polyhedron_faces:
                faces_id_list.InsertNextId(len(face))  # Number of vertices in this face
                for idx in face:
                    faces_id_list.InsertNextId(idx)

            grid.InsertNextCell(vtk.VTK_POLYHEDRON, faces_id_list)

        grid.SetPoints(points_vtk)

        if self.PostVelocity:
            velocity_array = vtk.vtkDoubleArray()
            velocity_array.SetName("Velocity")
            velocity_array.SetNumberOfComponents(3)

            i = 0
            for ii in range(len(self.polygon_centers_p)):
                for jj in range(len(self.polygon_origins_p[i])):
                    velocity_array.InsertNextTuple(self.polygon_velocity_p[i])
                i += 1

            grid.GetPointData().AddArray(velocity_array)

        if self.PostAngularVelocity:
            angular_velocity_array = vtk.vtkDoubleArray()
            angular_velocity_array.SetName("Angular Velocity")
            angular_velocity_array.SetNumberOfComponents(3)

            i = 0
            for ii in range(len(self.polygon_centers_p)):
                for jj in range(len(self.polygon_origins_p[i])):
                    angular_velocity_array.InsertNextTuple(self.polygon_angular_velocity_p[i])
                i += 1

            grid.GetPointData().AddArray(angular_velocity_array)

        if self.PostDisplacement:
            displacement_array = vtk.vtkDoubleArray()
            displacement_array.SetName("Displacement")
            displacement_array.SetNumberOfComponents(3)

            i = 0
            for ii in range(len(self.polygon_centers_p)):
                for jj in range(len(self.polygon_origins_p[i])):
                    displacement_array.InsertNextTuple(self.polygon_displacement_p[i])
                i += 1

            grid.GetPointData().AddArray(displacement_array)

        if self.PostTotalForces:
            totalforces_array = vtk.vtkDoubleArray()
            totalforces_array.SetName("Total Forces")
            totalforces_array.SetNumberOfComponents(3)

            i = 0
            for ii in range(len(self.polygon_centers_p)):
                for jj in range(len(self.polygon_origins_p[i])):
                    totalforces_array.InsertNextTuple(self.polygon_totalforces_p[i])
                i += 1

            grid.GetPointData().AddArray(totalforces_array)

        polygon_id_array = vtk.vtkIntArray()
        polygon_id_array.SetName("Polygon Id")
        polygon_id_array.SetNumberOfComponents(1)

        i = 0
        for ii in range(len(self.polygon_centers_p)):
            for jj in range(len(self.polygon_origins_p[i])):
                polygon_id_array.InsertNextTuple(self.polygon_id_p[i])
            i += 1

        grid.GetPointData().AddArray(polygon_id_array)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        file_path = os.path.join(self.vtk_post_path_directory, "polyhedron_p_{}.vtu".format(poly_output_cnt))
        writer.SetFileName(file_path)
        writer.SetInputData(grid)
        writer.Write()

    def WriteOutWalls(self, poly_output_cnt):

        polyhedron_points = []
        i = 0
        for polygon_center in self.polygon_centers_w:
            polygon = polygon_center + self.polygon_origins_w[i]
            for poly in polygon:
                polyhedron_points.append(poly)
            i += 1

        thick_wall = False
        for polyhedron_face in self.polyhedron_faces_w:
            if len(polyhedron_face) == 1:
                thick_wall = True
                break

        points_vtk = vtk.vtkPoints()
        cells_vtk = vtk.vtkCellArray()
        grid = vtk.vtkUnstructuredGrid()

        if thick_wall:
            polyhedron_faces = []
            i = 0
            offset = 0
            for polyhedron_face in self.polyhedron_faces_w:
                if i > 0:
                    offset += len(self.polygon_origins_w[i-1])
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

        else:
            polyhedron_faces_list = []
            polyhedron_faces = []
            i = 0
            offset = 0
            for polyhedron_face in self.polyhedron_faces_w:
                if i > 0:
                    offset += len(self.polygon_origins_w[i-1])
                for face in polyhedron_face:
                    adjusted_face = [v + offset for v in face]
                    polyhedron_faces.append(adjusted_face)
                polyhedron_faces_list.append(polyhedron_faces)
                polyhedron_faces = []
                i += 1
            
            for point_coord in polyhedron_points:
                points_vtk.InsertNextPoint(point_coord)

            for poly_id, polyhedron_faces in enumerate(polyhedron_faces_list):

                faces_id_list = vtk.vtkIdList()
                faces_id_list.InsertNextId(len(polyhedron_faces))  # Number of faces
                for face in polyhedron_faces:
                    faces_id_list.InsertNextId(len(face))  # Number of vertices in this face
                    for idx in face:
                        faces_id_list.InsertNextId(idx)

                grid.InsertNextCell(vtk.VTK_POLYHEDRON, faces_id_list)

            grid.SetPoints(points_vtk)

        if self.PostVelocity:
            velocity_array = vtk.vtkDoubleArray()
            velocity_array.SetName("Velocity")
            velocity_array.SetNumberOfComponents(3)

            i = 0
            for ii in range(len(self.polygon_centers_w)):
                for jj in range(len(self.polygon_origins_w[i])):
                    velocity_array.InsertNextTuple(self.polygon_velocity_w[i])
                i += 1

            grid.GetPointData().AddArray(velocity_array)

        if self.PostAngularVelocity:
            angular_velocity_array = vtk.vtkDoubleArray()
            angular_velocity_array.SetName("Angular Velocity")
            angular_velocity_array.SetNumberOfComponents(3)

            i = 0
            for ii in range(len(self.polygon_centers_w)):
                for jj in range(len(self.polygon_origins_w[i])):
                    angular_velocity_array.InsertNextTuple(self.polygon_angular_velocity_w[i])
                i += 1

            grid.GetPointData().AddArray(angular_velocity_array)

        if self.PostDisplacement:
            displacement_array = vtk.vtkDoubleArray()
            displacement_array.SetName("Displacement")
            displacement_array.SetNumberOfComponents(3)

            i = 0
            for ii in range(len(self.polygon_centers_w)):
                for jj in range(len(self.polygon_origins_w[i])):
                    displacement_array.InsertNextTuple(self.polygon_displacement_w[i])
                i += 1

            grid.GetPointData().AddArray(displacement_array)

        if self.PostTotalForces:
            totalforces_array = vtk.vtkDoubleArray()
            totalforces_array.SetName("Total Forces")
            totalforces_array.SetNumberOfComponents(3)

            i = 0
            for ii in range(len(self.polygon_centers_w)):
                for jj in range(len(self.polygon_origins_w[i])):
                    totalforces_array.InsertNextTuple(self.polygon_totalforces_w[i])
                i += 1

            grid.GetPointData().AddArray(totalforces_array)

        polygon_id_array = vtk.vtkIntArray()
        polygon_id_array.SetName("Polygon Id")
        polygon_id_array.SetNumberOfComponents(1)

        i = 0
        for ii in range(len(self.polygon_centers_w)):
            for jj in range(len(self.polygon_origins_w[i])):
                polygon_id_array.InsertNextTuple(self.polygon_id_w[i])
            i += 1

        grid.GetPointData().AddArray(polygon_id_array)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        file_path = os.path.join(self.vtk_post_path_directory, "polyhedron_w_{}.vtu".format(poly_output_cnt))
        writer.SetFileName(file_path)
        writer.SetInputData(grid)
        writer.Write()
    
    def WriteResults(self, poly_output_cnt):

        self.DataPreparation()

        if self.find_smp_p and self.find_smp_w:
            self.WriteOutParticles(poly_output_cnt)
            self.WriteOutWalls(poly_output_cnt)
        else:
            self.WriteOutAll(poly_output_cnt)
        
        
    