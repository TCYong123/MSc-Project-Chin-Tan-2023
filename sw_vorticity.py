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
Q = []

dmax = args.dmax

dt = [2.0**(-n) for n in range(-2,5)]
for i, T in enumerate(dt):
    if args.w:
        mesh = sw_create_mesh.main(["--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)])
        sw_im.main(["--write=2", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
        sw_trbdf2.main(["--write=2", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], mesh)

    with fd.CheckpointFile("vorticity_dt"+str(T*60*60)+"_"+str(dmax)+".h5", 'r') as afile:
        mesh = afile.load_mesh("Mesh")
        qI = afile.load_function(mesh, "q_outI")
        qT = afile.load_function(mesh, "q_outT")

    q_error = (fd.sqrt(fd.assemble(fd.dot(qT - qI, qT - qI) * fd.dx)))


    Q.append(q_error)
    print(f"L2 error at dt: {T} is vorticity: {Q[i]}")


plt.figure()
plt.plot(dt, Q, label='vorticity')
plt.xscale('log', base=2)
plt.yscale('log', base=2)
plt.legend()
plt.xlabel('dt')
plt.savefig('vorticity_'+str(dmax)+'.png')
        