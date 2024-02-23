from detran import Mesh1D

def get_core_mesh(cells_per_assembly=10) :
    """ Get the 5-assembly (with reflectors) mesh.  Check the cells per assembly
        is sane for the coarse mesh factor.
    """
    cm = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]
    fm = [cells_per_assembly]*7
    mt = [0, 1, 2, 3, 4, 5, 0]
    mesh = Mesh1D(fm, cm, mt)
    return mesh

def get_assembly_mesh(cells_per_assembly=10) :
    """ Get mesh for just one assembly. 
    """
    cm = [0.0, 10.0]
    fm = [cells_per_assembly]
    mt = [2]
    mesh = Mesh1D(fm, cm, mt)
    return mesh


if __name__ == "__main__":

    get_core_mesh(3).display