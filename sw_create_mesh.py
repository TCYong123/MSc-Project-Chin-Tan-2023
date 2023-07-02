import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import firedrake as fd
#get command arguments
from petsc4py import PETSc
PETSc.Sys.popErrorHandler()
import mg
import argparse

def main(raw_args=None):
    parser = argparse.ArgumentParser(description='Williamson 5 testcase for augmented Lagrangian solver.')
    parser.add_argument('--base_level', type=int, default=1, help='Base refinement level of icosahedral grid for MG solve. Default 1.')
    parser.add_argument('--ref_level', type=int, default=5, help='Refinement level of icosahedral grid. Default 5.')
    parser.add_argument('--dmax', type=float, default=15, help='Final time in days. Default 15.')
    parser.add_argument('--dumpt', type=float, default=24, help='Dump time in hours. Default 24.')
    parser.add_argument('--gamma', type=float, default=1.0e5, help='Augmented Lagrangian scaling parameter. Default 10000 for AL mode.')
    parser.add_argument('--solver_mode', type=str, default='monolithic', help='Solver strategy. monolithic=use monolithic MG with Schwarz smoothers. AL=use augmented Lagrangian formulation. Default = monolithic')
    parser.add_argument('--dt', type=float, default=1, help='Timestep in hours. Default 1.')
    parser.add_argument('--filename', type=str, default='w5aug')
    parser.add_argument('--coords_degree', type=int, default=1, help='Degree of polynomials for sphere mesh approximation.')
    parser.add_argument('--degree', type=int, default=1, help='Degree of finite element space (the DG space).')
    parser.add_argument('--kspschur', type=int, default=40, help='Max number of KSP iterations on the Schur complement. Default 40.')
    parser.add_argument('--kspmg', type=int, default=3, help='Max number of KSP iterations in the MG levels. Default 3.')
    parser.add_argument('--tlblock', type=str, default='mg', help='Solver for the velocity-velocity block. mg==Multigrid with patchPC, lu==direct solver with MUMPS, patch==just do a patch smoother. Default is mg')
    parser.add_argument('--schurpc', type=str, default='mass', help='Preconditioner for the Schur complement. mass==mass inverse, helmholtz==helmholtz inverse * laplace * mass inverse. Default is mass')
    parser.add_argument('--show_args', action='store_true', help='Output all the arguments.')
    parser.add_argument('--time_scheme', type=int, default=1, help='Timestepping scheme. 0=Crank-Nicholson. 1=Implicit midpoint rule.')

    args = parser.parse_known_args(raw_args)
    args = args[0]

    if args.show_args:
        PETSc.Sys.Print(args)

    # some domain, parameters and FS setup
    R0 = 6371220.
    H = fd.Constant(5960.)
    base_level = args.base_level
    nrefs = args.ref_level - base_level
    name = args.filename
    deg = args.coords_degree
    distribution_parameters = {"partition": True, "overlap_type": (fd.DistributedMeshOverlapType.VERTEX, 2)}
    #distribution_parameters = {"partition": True, "overlap_type": (fd.DistributedMeshOverlapType.FACET, 2)}


    def high_order_mesh_hierarchy(mh, degree, R0):
        meshes = []
        for m in mh:
            X = fd.VectorFunctionSpace(m, "Lagrange", degree)
            new_coords = fd.interpolate(m.coordinates, X)
            x, y, z = new_coords
            r = (x**2 + y**2 + z**2)**0.5
            new_coords.interpolate
            
            
            (R0*new_coords/r)
            new_mesh = fd.Mesh(new_coords)
            meshes.append(new_mesh)

        return fd.HierarchyBase(meshes, mh.coarse_to_fine_cells,
                                mh.fine_to_coarse_cells,
                                mh.refinements_per_level, mh.nested)

    if args.tlblock == "mg":
        basemesh = fd.IcosahedralSphereMesh(radius=R0,
                                            refinement_level=base_level,
                                            degree=1,
                                            distribution_parameters = distribution_parameters)
        del basemesh._radius
        mh = fd.MeshHierarchy(basemesh, nrefs)
        mh = high_order_mesh_hierarchy(mh, deg, R0)
        for mesh in mh:
            xf = mesh.coordinates
            mesh.transfer_coordinates = fd.Function(xf)
            x = fd.SpatialCoordinate(mesh)
            r = (x[0]**2 + x[1]**2 + x[2]**2)**0.5
            xf.interpolate(R0*xf/r)
            mesh.init_cell_orientations(x)
        mesh = mh[-1]
    else:
        mesh = fd.IcosahedralSphereMesh(radius=R0,
                                        refinement_level=args.ref_level, degree=deg,
                                        distribution_parameters = distribution_parameters)
        x = fd.SpatialCoordinate(mesh)
        mesh.init_cell_orientations(x)

    mesh.name = "Mesh"

    return mesh

if __name__ == "__main__":
    main()        