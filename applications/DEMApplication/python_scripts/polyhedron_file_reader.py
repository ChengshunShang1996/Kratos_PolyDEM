import KratosMultiphysics as Kratos


def ReadNextLine(f):
    while True:
        nextline=next(f)
        if nextline.startswith("//") == False :
            return nextline.split()

def ReadPolyhedronFile(filename):

    f = open(filename, 'r')
    list_of_vertices = []
    list_of_faces = []
    volume = []
    size = []

    for line in f:
        if line.startswith("//"):
            continue
        if line.startswith("Name"):
            data = ReadNextLine(f)
            name = data[0]
        if line.startswith("Begin Vertices"):
            while True:
                nextline=next(f)
                if nextline.startswith("//"):
                    continue
                if nextline.startswith("End Vertices"):
                    break
                vertices = nextline.strip().split()
                for vertex in vertices:
                    vertices = list(map(float, vertex.strip('()').split(',')))
                    list_of_vertices.append(vertices)
        if line.startswith("Begin Faces"):
            while True:
                nextline=next(f)
                if nextline.startswith("//"):
                    continue
                if nextline.startswith("End Faces"):
                    break
                faces = nextline.strip().split()
                for face in faces:
                    vertices = list(map(int, face.strip('()').split(',')))
                    list_of_faces.append(vertices)
        if line.startswith("Size"):
            data = ReadNextLine(f)
            size = [float(data[0])]
        if line.startswith("Volume"):
            data = ReadNextLine(f)
            volume = [float(data[0])]
    f.close()

    if len(volume)==0 or len(size)==0 or len(list_of_vertices)==0 or len(list_of_faces)==0 :
        message = "\n\n" + "************  ERROR!   Problems reading polyhedron file: " + filename + "   ***************\n\n"
        Kratos.Logger.PrintInfo(message)
    else:
        Kratos.Logger.PrintInfo("Polyhedron file "+ filename + " was read correctly")

    return [name, list_of_vertices, list_of_faces, size[0], volume[0]]