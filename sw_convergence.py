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

dt = [2.0**(-n) for n in range(-2, 10)]
for i, T in enumerate(dt):
    if args.w:
        mesh = sw_create_mesh.main(["--ref_level=3", "--dmax=1", "--dt="+str(T)])
        sw_im.main(["--write=1", "--ref_level=3", "--dmax=1", "--dt="+str(T)], mesh)
        sw_trbdf2.main(["--write=1", "--ref_level=3", "--dmax=1", "--dt="+str(T)], mesh)

    with fd.CheckpointFile("convergence_dt"+str(60*60*T)+".h5", 'r') as afile:
        mesh = afile.load_mesh("Mesh")
        uI = afile.load_function(mesh, "u_outI")
        hI = afile.load_function(mesh, "h_outI")
        uT = afile.load_function(mesh, "u_outT")
        hT = afile.load_function(mesh, "h_outT")

    u_error = (fd.sqrt(fd.assemble(fd.dot(uT - uI, uT - uI) * fd.dx)))
    h_error = (fd.sqrt(fd.assemble(fd.dot(hT - hI, hT - hI) * fd.dx)))

    U.append(u_error)
    H.append(h_error)
    print(f"L2 error at dt: {T} is u: {U[i]}, h: {H[i]}")

plt.figure()
plt.plot(dt, U, label='velocity')
plt.plot(dt, H, label='height')
plt.xscale('log', base=2)
plt.yscale('log', base=2)
plt.legend()
plt.xlabel('dt')
plt.savefig('convergence.png')
        