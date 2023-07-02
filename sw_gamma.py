import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import firedrake as fd
import matplotlib.pyplot as plt
import sw_im
import sw_trbdf2
import sw_create_mesh
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--w', type=bool, default=False) #Toggle writing of files
args = parser.parse_known_args()
args = args[0]
U = []
H = []

dmax = 1.0
a = 1
b = 9+1
gamma = ([round(0.1*n, 2)  for n in range(a, b)])
gamma.append(-1.0)


if args.w:
    mesh = sw_create_mesh.main(["--ref_level=5", "--dmax="+str(dmax), "--dt=1"])

    for i, G in enumerate(gamma):   
        sw_trbdf2.main(["--write_height=1", "--ref_level=5", "--dmax="+str(dmax), "--dt=1"], mesh, G)

    sw_im.main(["--write_height=1", "--ref_level=5", "--dmax="+str(dmax), "--dt=1"], mesh, gamma)

for i, G in enumerate(gamma):
    with fd.CheckpointFile("gamma_dt"+str(G)+"_"+str(dmax)+".h5", 'r') as afile:
        mesh = afile.load_mesh("Mesh")
        uI = afile.load_function(mesh, "u_outI")
        hI = afile.load_function(mesh, "h_outI")
        uT = afile.load_function(mesh, "u_outT")
        hT = afile.load_function(mesh, "h_outT")

    u_error = (fd.sqrt(fd.assemble(fd.dot(uT - uI, uT - uI) * fd.dx)))
    h_error = (fd.sqrt(fd.assemble(fd.dot(hT - hI, hT - hI) * fd.dx)))

    U.append(u_error)
    H.append(h_error)

    if G == -1.0:
        G = 2 - 2**0.5
    print(f"L2 error at gamma: {G} is u: {U[i]}, h: {H[i]}")

gamma.remove(-1.0)
gamma.append(2-2**0.5)

plt.figure()
plt.plot(gamma, U, 'x',label='velocity')
plt.plot(gamma, H, 'x',label='height')
plt.legend()
plt.xlabel('gamma')
plt.title('L2 error of IM and TRBDF2 for varying Gamma')
plt.savefig('gamma'+str(dmax)+'.png')
        