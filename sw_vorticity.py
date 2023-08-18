import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import firedrake as fd
import matplotlib.pyplot as plt
import sw_im
import sw_trbdf2
import sw_im_R
import sw_trbdf2_R
import sw_im_E
import sw_trbdf2_E
import sw_create_mesh
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--w', type=bool, default=False) #Toggle writing of files
parser.add_argument('--dmax', type=float, default=50)
parser.add_argument('--dt', type=float, default=0)
args = parser.parse_known_args()
args = args[0]
dmax = args.dmax
dt = 2.0**(-args.dt)

if args.w == True:
    mesh = sw_create_mesh.main(["--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)])
    sw_im.main(["--vorticity=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_trbdf2.main(["--vorticity=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_im_R.main(["--vorticity=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_trbdf2_R.main(["--vorticity=1",  "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_im_E.main(["--vorticity=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_trbdf2_E.main(["--vorticity=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)

T = dt*60*60
im_vorticity = np.loadtxt("im_vorticity"+str(T)+"_"+str(dmax)+".array")
tr_vorticity = np.loadtxt("tr_vorticity"+str(T)+"_"+str(dmax)+".array")
imR_vorticity = np.loadtxt("imR_vorticity"+str(T)+"_"+str(dmax)+".array")
trR_vorticity= np.loadtxt("trR_vorticity"+str(T)+"_"+str(dmax)+".array")
imE_vorticity = np.loadtxt("imE_vorticity"+str(T)+"_"+str(dmax)+".array")
trE_vorticity= np.loadtxt("trE_vorticity"+str(T)+"_"+str(dmax)+".array")

n = np.arange(len(im_vorticity)*dt, step=dt)

plt.figure()
plt.semilogy(n, im_vorticity, label="IM")
plt.semilogy(n, tr_vorticity, label="TR-BDF2")
plt.semilogy(n, imR_vorticity, label="IM(R)")
plt.semilogy(n, trR_vorticity, label="TR-BDF2(R)")
plt.semilogy(n, imE_vorticity, label="IM(E)")
plt.semilogy(n, trE_vorticity, label="TR-BDF2(E)")
plt.legend()
plt.xlabel("Time (hours)")
plt.ylabel("Vorticity level")
plt.title("Vorticity comparison for different methods, dt="+str(dt))
plt.savefig('vorticity_'+str(T)+'_'+str(dmax)+'.png')

#Potential Vorticity:  1/h0 * (fd.dot(perp(fd.grad), u0) + f)


# import os
# os.environ["OMP_NUM_THREADS"] = "1"
# os.environ["OPENBLAS_NUM_THREADS"] = "1" 
# import firedrake as fd
# import matplotlib.pyplot as plt
# import sw_im
# import sw_trbdf2
# import sw_create_mesh
# import argparse
# import numpy as np

# parser = argparse.ArgumentParser()
# parser.add_argument('--w', type=bool, default=False) #Toggle writing of files
# parser.add_argument('--dmax', type=float, default=15.0)
# args = parser.parse_known_args()
# args = args[0]
# Q = []

# dmax = args.dmax

# dt = [2.0**(-n) for n in range(-2,5)]
# for i, T in enumerate(dt):
#     if args.w:
#         mesh = sw_create_mesh.main(["--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)])
#         sw_im.main(["--write=2", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
#         sw_trbdf2.main(["--write=2", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], mesh)

#     with fd.CheckpointFile("vorticity_dt"+str(T*60*60)+"_"+str(dmax)+".h5", 'r') as afile:
#         mesh = afile.load_mesh("Mesh")
#         qI = afile.load_function(mesh, "q_outI")
#         qT = afile.load_function(mesh, "q_outT")

#     q_error = (fd.sqrt(fd.assemble(fd.dot(qT - qI, qT - qI) * fd.dx)))


#     Q.append(q_error)
#     print(f"L2 error at dt: {T} is vorticity: {Q[i]}")


# plt.figure()
# plt.plot(dt, Q, label='vorticity')
# plt.xscale('log', base=2)
# plt.yscale('log', base=2)
# plt.legend()
# plt.xlabel('dt')
# plt.savefig('vorticity_'+str(dmax)+'.png')
        