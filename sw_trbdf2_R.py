import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import firedrake as fd
#get command arguments
from petsc4py import PETSc
PETSc.Sys.popErrorHandler()
import mg
import argparse
import numpy as np

def main(raw_args=None, mesh=None, Gam=None):
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
    parser.add_argument('--time_scheme', type=int, default=0, help='Timestepping scheme. 0=TR-BDF2.')
    parser.add_argument('--write', type=int, default=0, help='Write files for convergence (dt varies). 0=None, 1=Convergence(IM), 2=Vorticity, 3=Rosenbrock')
    parser.add_argument('--array', type=int, default=0, help='Write array. 0=None, 1=height, 2=vorticity')

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

    if mesh == None:
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

    R0 = fd.Constant(R0)
    cx, cy, cz = fd.SpatialCoordinate(mesh)

    outward_normals = fd.CellNormal(mesh)


    def perp(u):
        return fd.cross(outward_normals, u)


    degree = args.degree
    V1 = fd.FunctionSpace(mesh, "BDM", degree+1)
    V2 = fd.FunctionSpace(mesh, "DG", degree)
    V0 = fd.FunctionSpace(mesh, "CG", degree+2)
    W = fd.MixedFunctionSpace((V1, V2))

    u, eta = fd.TrialFunctions(W)
    v, phi = fd.TestFunctions(W)

    Omega = fd.Constant(7.292e-5)  # rotation rate
    f = 2*Omega*cz/fd.Constant(R0)  # Coriolis parameter
    g = fd.Constant(9.8)  # Gravitational constant
    b = fd.Function(V2, name="Topography")
    c = fd.sqrt(g*H)
    if args.solver_mode == "AL":
        gamma0 = args.gamma
        gamma = fd.Constant(gamma0)
    else:
        gamma0 = 0.
        gamma = fd.Constant(gamma0)

    # D = eta + b

    One = fd.Function(V2).assign(1.0)

    u, eta = fd.TrialFunctions(W)
    v, phi = fd.TestFunctions(W)

    dx = fd.dx

    Un = fd.Function(W, name="Un")
    Ug = fd.Function(W, name="Ug")
    Unp1 = fd.Function(W, name="Unp1")


    u0, h0 = fd.split(Un)
    ug, hg = fd.split(Ug)
    u1, h1 = fd.split(Unp1)
    n = fd.FacetNormal(mesh)


    def both(u):
        return 2*fd.avg(u)


    dT = fd.Constant(0.)
    dS = fd.dS


    def u_op(v, u, h):
        Upwind = 0.5 * (fd.sign(fd.dot(u, n)) + 1)
        K = 0.5*fd.inner(u, u)
        return (fd.inner(v, f*perp(u))*dx
                - fd.inner(perp(fd.grad(fd.inner(v, perp(u)))), u)*dx
                + fd.inner(both(perp(n)*fd.inner(v, perp(u))),
                            both(Upwind*u))*dS
                - fd.div(v)*(g*(h + b) + K)*dx)


    def h_op(phi, u, h):
        uup = 0.5 * (fd.dot(u, n) + abs(fd.dot(u, n)))
        return (- fd.inner(fd.grad(phi), u)*h*dx
                + fd.jump(phi)*(uup('+')*h('+')
                                - uup('-')*h('-'))*dS)



    if args.time_scheme == 0:
        "TR-BDF2"
        if Gam == None:
            gam = fd.Constant(2-(2**0.5))
        elif Gam == -1:
            gam = fd.Constant(2-(2**0.5))
        else: 
            gam = Gam
        half = fd.Constant(0.5)
        quo1 = fd.Constant(1/(gam*(2-gam)))
        quo2 = fd.Constant(1/(2-gam))
        con1 = fd.Constant((1-gam)**2*quo1)
        con2 = fd.Constant((1-gam)*quo2)

        testeqn1 = (
            fd.inner(v, ug - u0)*dx
            + gam*half*dT*u_op(v, u0, h0)
            + gam*half*dT*u_op(v, ug, hg)
            + phi*(hg - h0)*dx
            + gam*half*dT*h_op(phi, u0, h0)
            + gam*half*dT*h_op(phi, ug, hg))

        eqn1 = testeqn1 \
            + 0    

        testeqn2 = (
            fd.inner(v, u1 - quo1*ug + con1*u0)*dx
            + con2*dT*u_op(v, u1, h1)
            + phi*(h1 - quo1*hg + con1*h0)*dx
            + con2*dT*h_op(phi, u1, h1))         
        
        eqn2 = testeqn2 \
            + 0  


    else:
        raise NotImplementedError


    # TRAPEZOIDAL RULE
    # U^{n+1} - U^n + dt*( N(U^{n+1}) + N(U^n) )/2 = 0.

    # BDF2
    # y^{n+2} - 4/3*y^{n+1} + 1/3*y^n - 2/3*dt*( f(y^{n+2})) = 0 

    # TR-BDF2
    # y^{n+g} - y^n - g*dt/2*( f(y^{n+g}+f(y^{n}))) = 0   
    # y^{n+1} - 1/(g*(2-g))*y^{n+g} - (1-g)/(2-g)*dt*f(y^{n+1}) + (1-g)^2/(g*(2-g))*y^n = 0




    class HelmholtzPC(fd.PCBase):
        def initialize(self, pc):
            if pc.getType() != "python":
                raise ValueError("Expecting PC type python")
            prefix = pc.getOptionsPrefix() + "helmholtz_"

            mm_solve_parameters = {
                'ksp_type':'preonly',
                'pc_type':'bjacobi',
                'sub_pc_type':'lu',
            }

            # we assume P has things stuffed inside of it
            _, P = pc.getOperators()
            context = P.getPythonContext()
            appctx = context.appctx
            self.appctx = appctx

            # FunctionSpace checks
            u, v = context.a.arguments()
            if u.function_space() != v.function_space():
                raise ValueError("Pressure space test and trial space differ")

            # the mass solve
            a = u*v*fd.dx
            self.Msolver = fd.LinearSolver(fd.assemble(a),
                                        solver_parameters=
                                        mm_solve_parameters)
            # the Helmholtz solve
            eta0 = appctx.get("helmholtz_eta", 20)
            def get_laplace(q,phi):
                h = fd.avg(fd.CellVolume(mesh))/fd.FacetArea(mesh)
                mu = eta0/h
                n = fd.FacetNormal(mesh)
                ad = (- fd.jump(phi)*fd.jump(fd.inner(fd.grad(q),n))
                    - fd.jump(q)*fd.jump(fd.inner(fd.grad(phi),n)))*fd.dS
                ad +=  mu * fd.jump(phi)*fd.jump(q)*fd.dS
                ad += fd.inner(fd.grad(q), fd.grad(phi)) * fd.dx
                return ad

            a = (fd.Constant(2)/dT/H)*u*v*fd.dx + fd.Constant(0.5)*g*dT*get_laplace(u, v)
            #input and output functions
            V = u.function_space()
            self.xfstar = fd.Function(V) # the input residual
            self.xf = fd.Function(V) # the output function from Riesz map
            self.yf = fd.Function(V) # the preconditioned residual

            L = get_laplace(u, self.xf*gamma)
            hh_prob = fd.LinearVariationalProblem(a, L, self.yf)
            self.hh_solver = fd.LinearVariationalSolver(
                hh_prob,
                options_prefix=prefix)

        def update(self, pc):
            pass

        def applyTranspose(self, pc, x, y):
            raise NotImplementedError

        def apply(self, pc, x, y):

            # copy petsc vec into Function
            with self.xfstar.dat.vec_wo as v:
                x.copy(v)
            
            #do the mass solver, solve(x, b)
            self.Msolver.solve(self.xf, self.xfstar)

            # get the mean
            xbar = fd.assemble(self.xf*fd.dx)/fd.assemble(One*fd.dx)
            self.xf -= xbar
            
            #do the Helmholtz solver
            self.hh_solver.solve()

            # add the mean
            self.yf += xbar/2*dT*gamma*H
            
            # copy petsc vec into Function
            with self.yf.dat.vec_ro as v:
                v.copy(y)


    if args.solver_mode == 'AL':
        
        sparameters = {
            "mat_type":"matfree",
            'snes_monitor': None,
            "ksp_type": "fgmres",
            "ksp_gmres_modifiedgramschmidt": None,
            'ksp_monitor': None,
            "pc_type": "fieldsplit",
            "pc_fieldsplit_type": "schur",
            #"pc_fieldsplit_schur_fact_type": "full",
            "pc_fieldsplit_off_diag_use_amat": True,
        }

        bottomright_helm = {
            "ksp_type": "fgmres",
            "ksp_monitor": None,
            "ksp_gmres_modifiedgramschmidt": None,
            "ksp_max_it": args.kspschur,
            "pc_type": "python",
            "pc_python_type": "__main__.HelmholtzPC",
            "helmholtz" :
            {"ksp_type":"preonly",
            "pc_type": "lu"}
        }

        bottomright_mass = {
            "ksp_type": "preonly",
            #"ksp_monitor":None,
            "ksp_gmres_modifiedgramschmidt": None,
            "ksp_max_it": args.kspschur,
            #"ksp_monitor":None,
            "pc_type": "python",
            "pc_python_type": "firedrake.MassInvPC",
            "Mp_pc_type": "bjacobi",
            "Mp_sub_pc_type": "ilu"
        }

        if args.schurpc == "mass":
            sparameters["fieldsplit_1"] = bottomright_mass
        elif args.schurpc == "helmholtz":
            sparameters["fieldsplit_1"] = bottomright_helm
        else:
            raise KeyError('Unknown Schur PC option.')

        topleft_LU = {
            "ksp_type": "preonly",
            "pc_type": "python",
            "pc_python_type": "firedrake.AssembledPC",
            "assembled_pc_type": "lu",
            "assembled_pc_factor_mat_solver_type": "mumps"
        }

        topleft_MG = {
            "ksp_type": "preonly",
            "pc_type": "mg",
            #"pc_mg_type": "full",
            "mg_coarse_ksp_type": "preonly",
            "mg_coarse_pc_type": "python",
            "mg_coarse_pc_python_type": "firedrake.AssembledPC",
            "mg_coarse_assembled_pc_type": "lu",
            "mg_coarse_assembled_pc_factor_mat_solver_type": "mumps",
            "mg_levels_ksp_type": "gmres",
            "mg_levels_ksp_max_it": args.kspmg,
            "mg_levels_pc_type": "python",
            "mg_levels_pc_python_type": "firedrake.PatchPC",
            "mg_levels_patch_pc_patch_save_operators": True,
            "mg_levels_patch_pc_patch_partition_of_unity": False,
            "mg_levels_patch_pc_patch_sub_mat_type": "seqaij",
            "mg_levels_patch_pc_patch_construct_type": "star",
            "mg_levels_patch_pc_patch_multiplicative": False,
            "mg_levels_patch_pc_patch_symmetrise_sweep": False,
            "mg_levels_patch_pc_patch_construct_dim": 0,
            "mg_levels_patch_pc_patch_sub_mat_type": "seqdense",
            "mg_levels_patch_pc_patch_dense_inverse": True,
            "mg_levels_patch_pc_patch_precompute_element_tensors": True,
            "mg_levels_patch_sub_pc_factor_mat_solver_type": "petsc",
            "mg_levels_patch_sub_ksp_type": "preonly",
            "mg_levels_patch_sub_pc_type": "lu",
        }

        topleft_MGs = {
            "ksp_type": "preonly",
            "ksp_max_it": 3,
            "pc_type": "mg",
            "mg_coarse_ksp_type": "preonly",
            "mg_coarse_pc_type": "python",
            "mg_coarse_pc_python_type": "firedrake.AssembledPC",
            "mg_coarse_assembled_pc_type": "lu",
            "mg_coarse_assembled_pc_factor_mat_solver_type": "mumps",
            "mg_levels_ksp_type": "gmres",
            "mg_levels_ksp_max_it": args.kspmg,
            "mg_levels_pc_type": "python",
            "mg_levels_pc_python_type": "firedrake.AssembledPC",
            "mg_levels_assembled_pc_type": "python",
            "mg_levels_assembled_pc_python_type": "firedrake.ASMStarPC",
            "mg_levels_assembled_pc_star_backend": "tinyasm",
            "mg_levels_assmbled_pc_star_construct_dim": 0
        }
        
        topleft_smoother = {
            "ksp_type": "gmres",
            "ksp_max_it": 3,
            "ksp_monitor": None,
            "pc_type": "python",
            "pc_python_type": "firedrake.PatchPC",
            "patch_pc_patch_save_operators": True,
            "patch_pc_patch_partition_of_unity": False,
            "patch_pc_patch_sub_mat_type": "seqaij",
            "patch_pc_patch_construct_type": "star",
            "patch_pc_patch_multiplicative": False,
            "patch_pc_patch_symmetrise_sweep": False,
            "patch_pc_patch_construct_dim": 0,
            "patch_sub_ksp_type": "preonly",
            "patch_sub_pc_type": "lu",
        }

        if args.tlblock == "mg":
            sparameters["fieldsplit_0"] = topleft_MG
        elif args.tlblock == "patch":
            sparameters["fieldsplit_0"] = topleft_smoother
        else:
            assert(args.tlblock=="lu")
            sparameters["fieldsplit_0"] = topleft_LU
    elif args.solver_mode == 'monolithic':
        # monolithic solver options

        sparameters = {
            "snes_monitor": None,
            "snes_type":"ksponly",
            "mat_type": "matfree",
            "ksp_type": "fgmres",
            "ksp_monitor_true_residual": None,
            "ksp_converged_reason": None,
            "ksp_atol": 1e-8,
            "ksp_rtol": 1e-8,
            "ksp_max_it": 400,
            "pc_type": "mg",
            "pc_mg_cycle_type": "v",
            "pc_mg_type": "multiplicative",
            "mg_levels_ksp_type": "gmres",
            "mg_levels_ksp_max_it": 3,
            #"mg_levels_ksp_convergence_test": "skip",
            "mg_levels_pc_type": "python",
            "mg_levels_pc_python_type": "firedrake.PatchPC",
            "mg_levels_patch_pc_patch_save_operators": True,
            "mg_levels_patch_pc_patch_partition_of_unity": True,
            "mg_levels_patch_pc_patch_sub_mat_type": "seqdense",
            "mg_levels_patch_pc_patch_construct_dim": 0,
            "mg_levels_patch_pc_patch_construct_type": "vanka",
            "mg_levels_patch_pc_patch_local_type": "additive",
            "mg_levels_patch_pc_patch_precompute_element_tensors": True,
            "mg_levels_patch_pc_patch_symmetrise_sweep": False,
            "mg_levels_patch_sub_ksp_type": "preonly",
            "mg_levels_patch_sub_pc_type": "lu",
            "mg_levels_patch_sub_pc_factor_shift_type": "nonzero",
            "mg_coarse_pc_type": "python",
            "mg_coarse_pc_python_type": "firedrake.AssembledPC",
            "mg_coarse_assembled_pc_type": "lu",
            "mg_coarse_assembled_pc_factor_mat_solver_type": "mumps",
        }
        
    dt = 60*60*args.dt
    dT.assign(dt)
    t = 0.

    nprob1 = fd.NonlinearVariationalProblem(eqn1, Ug)
    ctx = {"mu": gamma*2/g/dt}
    nsolver1 = fd.NonlinearVariationalSolver(nprob1,
                                            solver_parameters=sparameters,
                                            appctx=ctx)

    nprob2 = fd.NonlinearVariationalProblem(eqn2, Unp1)
    ctx = {"mu": gamma*2/g/dt}
    nsolver2 = fd.NonlinearVariationalSolver(nprob2,
                                            solver_parameters=sparameters,
                                            appctx=ctx)

    vtransfer = mg.ManifoldTransfer()
    tm = fd.TransferManager()
    transfers = {
        V1.ufl_element(): (vtransfer.prolong, vtransfer.restrict,
                        vtransfer.inject),
        V2.ufl_element(): (vtransfer.prolong, vtransfer.restrict,
                        vtransfer.inject)
    }
    transfermanager = fd.TransferManager(native_transfers=transfers)
    nsolver1.set_transfer_manager(transfermanager)
    nsolver2.set_transfer_manager(transfermanager)

    dmax = args.dmax
    hmax = 24*dmax
    tmax = 60.*60.*hmax
    hdump = args.dumpt
    dumpt = hdump*60.*60.
    tdump = 0.

    x = fd.SpatialCoordinate(mesh)
    u_0 = 20.0  # maximum amplitude of the zonal wind [m/s]
    u_max = fd.Constant(u_0)
    u_expr = fd.as_vector([-u_max*x[1]/R0, u_max*x[0]/R0, 0.0])
    eta_expr = - ((R0 * Omega * u_max + u_max*u_max/2.0)*(x[2]*x[2]/(R0*R0)))/g
    un = fd.Function(V1, name="Velocity").project(u_expr)
    etan = fd.Function(V2, name="Elevation").project(eta_expr)

    # Topography.
    rl = fd.pi/9.0
    lambda_x = fd.atan_2(x[1]/R0, x[0]/R0)
    lambda_c = -fd.pi/2.0
    phi_x = fd.asin(x[2]/R0)
    phi_c = fd.pi/6.0
    minarg = fd.min_value(pow(rl, 2),
                    pow(phi_x - phi_c, 2) + pow(lambda_x - lambda_c, 2))
    bexpr = 2000.0*(1 - fd.sqrt(minarg)/rl)
    b.interpolate(bexpr)

    u0, h0 = Un.subfunctions
    u0.assign(un)
    h0.assign(etan + H - b)

    q = fd.TrialFunction(V0)
    p = fd.TestFunction(V0)

    qn = fd.Function(V0, name="Relative Vorticity")
    veqn = q*p*dx + fd.inner(perp(fd.grad(p)), un)*dx
    vprob = fd.LinearVariationalProblem(fd.lhs(veqn), fd.rhs(veqn), qn)
    qparams = {'ksp_type':'cg'}
    qsolver = fd.LinearVariationalSolver(vprob,
                                        solver_parameters=qparams)

    file_sw = fd.File('tr'+name+'.pvd')
    etan.assign(h0 - H + b)
    un.assign(u0)
    qsolver.solve()
    file_sw.write(un, etan, qn)
    Unp1.assign(Un)

    V1DG = fd.VectorFunctionSpace(mesh,"DG",degree+1)
    u_out = fd.Function(V1DG, name="u_outR")
    h_out = fd.Function(V2, name="h_outR")
    q_out = fd.Function(V0, name="q_outR")

    ht_array = np.array([])
    vt_array = np.array([])
    qt_array = np.array([])
    b_array = np.array([])

    PETSc.Sys.Print('tmax', tmax, 'dt', dt)
    itcount = 0
    stepcount = 0
    while t < tmax + 0.5*dt:
        PETSc.Sys.Print(t)
        t += dt
        tdump += dt

        nsolver1.solve()
        nsolver2.solve()
        Un.assign(Unp1)
        res = fd.assemble(testeqn1+testeqn2)
        PETSc.Sys.Print(res.dat.data[0].max(), res.dat.data[0].min(),
            res.dat.data[1].max(), res.dat.data[1].min())          
        
        if tdump > dumpt - dt*0.5:
            etan.assign(h0 - H + b)
            un.assign(u0)
            qsolver.solve()
            file_sw.write(un, etan, qn)
            tdump -= dumpt
        stepcount += 1
        itcount += nsolver1.snes.getLinearSolveIterations()

        ht_array = np.append(ht_array, h0.dat.data[0])
        vt_array = np.append(vt_array, u0.dat.data[0])
        qt_array = np.append(qt_array, qn.dat.data[0])
        b_array = np.append(b_array, b.dat.data[0])


        if args.array == 1:
            if Gam == None:
                np.savetxt("trR_ht"+str(dt)+"_"+str(dmax)+".array", ht_array)
                np.savetxt("trR_vt"+str(dt)+"_"+str(dmax)+".array", vt_array)
                np.savetxt("trR_b"+str(dt)+"_"+str(dmax)+".array", b_array)
            else:
                np.savetxt("trR_ht"+str(dt)+"_"+str(dmax)+"_"+str(Gam)+".array", ht_array)
                np.savetxt("trR_vt"+str(dt)+"_"+str(dmax)+"_"+str(Gam)+".array", vt_array)
                np.savetxt("trR_b"+str(dt)+"_"+str(dmax)+".array", b_array)

        elif args.array == 2:
            np.savetxt("trR_vor"+str(dt)+"_"+str(dmax)+".array", qt_array) 

    PETSc.Sys.Print("Iterations", itcount, "its per step", itcount/stepcount,
                    "dt", dt, "tlblock", args.tlblock, "ref_level", args.ref_level, "dmax", args.dmax)
    

    # To be used with sw_im as file is appended after
    if args.write == 1:
        u_out.interpolate(u0)
        h_out.interpolate(h0)
        with fd.CheckpointFile("convergence_dt"+str(dt)+"_"+str(dmax)+".h5", 'a') as afile:
            afile.save_function(u_out)
            afile.save_function(h_out)

    elif args.write == 2:
        q_out.interpolate(qn)
        with fd.CheckpointFile("vorticity_dt"+str(dt)+"_"+str(dmax)+".h5", 'w') as afile:
            afile.save_function(q_out)

    elif args.write == 3:         
        u_out.interpolate(u0)
        h_out.interpolate(h0)
        with fd.CheckpointFile("rosenbrock_dt"+str(dt)+"_"+str(dmax)+".h5", 'w') as afile:
            afile.save_function(u_out)
            afile.save_function(h_out)


    if Gam != None:
        u_out.interpolate(u0)
        h_out.interpolate(h0)
        with fd.CheckpointFile("gamma_dt"+str(Gam)+"_"+str(dmax)+".h5", 'w') as afile:
            afile.save_function(u_out)
            afile.save_function(h_out)

if __name__ == "__main__":
    main()