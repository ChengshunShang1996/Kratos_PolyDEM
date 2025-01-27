#/////////////////////////////////////////////////////////////
#// Main author: Chengshun Shang (CIMNE)
#// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
#// Date: July 2024
#/////////////////////////////////////////////////////////////

import KratosMultiphysics as Kratos


def ReadNextLine(f):
    while True:
        nextline=next(f)
        if nextline.startswith("//") == False :
            return nextline.split()

def ReadPolyhedronFile(filename):

    f = open(filename, 'r')
    list_of_vertices_list = []
    list_of_faces_list = []
    list_of_size = []
    list_of_volume = []
    list_of_inertia_per_unit_mass_list = []

    for line in f:
        if line.startswith("//"):
            continue
        
        if line.startswith("Name"):
            data = ReadNextLine(f)
            name = data[0]
        
        if line.startswith("Begin Vertices"):
            while True:
                nextline=next(f)
                list_of_vertices = []
                if nextline.startswith("//"):
                    continue
                if nextline.startswith("End Vertices"):
                    break
                vertices_list = nextline.strip().split()
                for vertex in vertices_list:
                    vertices = list(map(float, vertex.strip('()').split(',')))
                    list_of_vertices.append(vertices)
                list_of_vertices_list.append(list_of_vertices)
        
        if line.startswith("Begin Faces"):
            while True:
                nextline=next(f)
                list_of_faces = []
                if nextline.startswith("//"):
                    continue
                if nextline.startswith("End Faces"):
                    break
                faces = nextline.strip().split()
                for face in faces:
                    vertices = list(map(int, face.strip('()').split(',')))
                    list_of_faces.append(vertices)
                list_of_faces_list.append(list_of_faces)
        
        if line.startswith("Begin Size"):
            while True:
                nextline=next(f)
                if nextline.startswith("//"):
                    continue
                if nextline.startswith("End Size"):
                    break
                list_of_size.append(float(nextline.split()[0]))

        if line.startswith("Begin Volume"):
            while True:
                nextline=next(f)
                if nextline.startswith("//"):
                    continue
                if nextline.startswith("End Volume"):
                    break
                list_of_volume.append(float(nextline.split()[0]))

        if line.startswith("Begin InertiaPerUnitMass"):
            while True:
                nextline=next(f)
                list_of_inertia_per_unit_mass = []
                if nextline.startswith("//"):
                    continue
                if nextline.startswith("End InertiaPerUnitMass"):
                    break
                inertia_per_unit_mass_list = nextline.split()[0]
                inertia = list(map(float, inertia_per_unit_mass_list.strip('()').split(',')))
                list_of_inertia_per_unit_mass.append(inertia)
                list_of_inertia_per_unit_mass_list.append(list_of_inertia_per_unit_mass)
    
    f.close()

    if len(list_of_volume)==0 or len(list_of_size)==0 or len(list_of_faces_list)==0 or len(list_of_vertices_list)==0 or len(list_of_inertia_per_unit_mass_list)==0:
        message = "\n\n" + "************  ERROR!   Problems reading polyhedron file: " + filename + "   ***************\n\n"
        Kratos.Logger.PrintInfo(message)
    else:
        Kratos.Logger.PrintInfo("Polyhedron file "+ filename + " was read correctly")

    return [name, list_of_vertices_list, list_of_faces_list, list_of_size, list_of_volume, list_of_inertia_per_unit_mass_list]