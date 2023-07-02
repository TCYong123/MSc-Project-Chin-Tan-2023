import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import firedrake as fd
import matplotlib.pyplot as plt
import sw_im
import sw_trbdf2
import sw_create_mesh
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--w', type=bool, default=False) #Toggle writing of files
parser.add_argument('--dmax', type=float, default=15.0)
args = parser.parse_known_args()
args = args[0]
U = []
H = []
dmax = args.dmax

dt = [2.0**(-n) for n in range(-2,5)]
for i, T in enumerate(dt):
    if args.w:
        mesh = sw_create_mesh.main(["--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)])
        sw_im.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
        sw_trbdf2.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], mesh)

    with fd.CheckpointFile("convergence_dt"+str(T*60*60)+"_"+str(dmax)+".h5", 'r') as afile:
        mesh = afile.load_mesh("Mesh")
        uI = afile.load_function(mesh, "u_outI")
        hI = afile.load_function(mesh, "h_outI")
        uT = afile.load_function(mesh, "u_outT")
        hT = afile.load_function(mesh, "h_outT")

    u_error = (fd.sqrt(fd.assemble(fd.dot(uT - uI, uT - uI) * fd.dx))) 
    h_error = (fd.sqrt(fd.assemble(fd.dot(hT - hI, hT - hI) * fd.dx)))

    if dt == 1.0:
        u_ref = u_error
        h_ref = h_error


    U.append(u_error)
    H.append(h_error)
    print(f"L2 error at dt: {T} is u: {U[i]}, h: {H[i]}")

U = np.array(U)
H = np.array(H)

plt.figure()
plt.plot(dt, U/u_error, label='normalized velocity')
plt.plot(dt, H/h_error, label='normalized height')
plt.xscale('log', base=2)
plt.yscale('log', base=2)
plt.legend()
plt.xlabel('dt')
plt.savefig('convergence_'+str(dmax)+'.png')
        