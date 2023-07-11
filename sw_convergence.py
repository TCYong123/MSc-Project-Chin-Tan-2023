import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import firedrake as fd
import matplotlib.pyplot as plt
import sw_im
import sw_trbdf2
import sw_im_R
import sw_trbdf2_R
import sw_create_mesh
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--w', type=bool, default=False) #Toggle writing of files
parser.add_argument('--dmax', type=float, default=15)
args = parser.parse_known_args()
args = args[0]
U = np.array([])
H = np.array([])
UI = np.array([])
HI = np.array([])
UR = np.array([])
HR = np.array([])
dmax = args.dmax

dt = [2.0**(-n) for n in range(0,12)]
for i, T in enumerate(dt):
    if args.w:
        mesh = sw_create_mesh.main(["--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)])
        sw_im.main(["--write=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
        sw_im_R.main(["--write=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
        sw_trbdf2_R.main(["--write=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
        sw_trbdf2.main(["--write=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], mesh)

    with fd.CheckpointFile("convergence_dt"+str(T*60*60)+"_"+str(dmax)+".h5", 'r') as afile:
        mesh = afile.load_mesh("Mesh")
        uI = afile.load_function(mesh, "u_outI")
        hI = afile.load_function(mesh, "h_outI")
        uT = afile.load_function(mesh, "u_outT")
        hT = afile.load_function(mesh, "h_outT")
        uIR = afile.load_function(mesh, "u_outIR")
        hIR = afile.load_function(mesh, "h_outIR")
        uTR = afile.load_function(mesh, "u_outTR")
        hTR = afile.load_function(mesh, "h_outTR")        

    u_error = (fd.sqrt(fd.assemble(fd.dot(uT - uI, uT - uI) * fd.dx))) 
    h_error = (fd.sqrt(fd.assemble(fd.dot(hT - hI, hT - hI) * fd.dx)))
    uI_error = (fd.sqrt(fd.assemble(fd.dot(uIR - uI, uIR - uI) * fd.dx))) 
    hI_error = (fd.sqrt(fd.assemble(fd.dot(hIR - hI, hIR - hI) * fd.dx)))
    uR_error = (fd.sqrt(fd.assemble(fd.dot(uTR - uI, uTR - uI) * fd.dx))) 
    hR_error = (fd.sqrt(fd.assemble(fd.dot(hTR - hI, hTR - hI) * fd.dx)))

    if dt == 1.0:
        u_ref = u_error
        h_ref = h_error
        uI_ref = uI_error
        hI_ref = hI_error
        uR_ref = uR_error
        hR_ref = hR_error

    U = np.append(U, u_error)
    H = np.append(H, h_error)
    UI = np.append(UI, uI_error)
    HI = np.append(HI, hI_error)
    UR = np.append(UR, uR_error)
    HR = np.append(HR, hR_error)
      
    print(f"L2 error at dt: {T} is u: {U[i]}, h: {H[i]}")

plt.figure()
plt.plot(dt, U/u_ref, label='normalized velocity TR-BDF2')
plt.plot(dt, H/h_ref, label='normalized height TR-BDF2')
plt.plot(dt, UI/uI_ref, label='normalized velocity IM(R)')
plt.plot(dt, HI/hI_ref, label='normalized height IM(R)')
plt.plot(dt, UR/uR_ref, label='normalized velocity TR-BDF2(R)')
plt.plot(dt, HR/hR_ref, label='normalized height TR-BDF2(R)')
# plt.xscale('log', base=2)
# plt.yscale('log', base=2)
plt.legend()
plt.xlabel('dt')
plt.title("Normalized error between methods and IM")
plt.savefig('convergence_'+str(dmax)+'.png')
        