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
parser.add_argument('--dmax', type=float, default=15.0)
args = parser.parse_known_args()
args = args[0]
U = np.array([])
H = np.array([])
UIR = np.array([])
HIR = np.array([])
UTR = np.array([])
HTR = np.array([])
UIE = np.array([])
HIE = np.array([])
UTE = np.array([])
HTE = np.array([])
dmax = args.dmax

dt = [2.0**(-n) for n in range(0,8)]
for i, T in enumerate(dt):
    if args.w:
        mesh = sw_create_mesh.main(["--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)])
        sw_im.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
        sw_trbdf2_R.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
        sw_im_R.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
        sw_trbdf2.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
        sw_im_E.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
        sw_trbdf2_E.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], mesh)

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
        uIE = afile.load_function(mesh, "u_outIE")
        hIE = afile.load_function(mesh, "h_outIE")
        uTE = afile.load_function(mesh, "u_outTE")
        hTE = afile.load_function(mesh, "h_outTE") 

    u_error = (fd.sqrt(fd.assemble(fd.dot(uT - uI, uT - uI) * fd.dx))) 
    h_error = (fd.sqrt(fd.assemble(fd.dot(hT - hI, hT - hI) * fd.dx)))
    uIR_error = (fd.sqrt(fd.assemble(fd.dot(uIR - uI, uIR - uI) * fd.dx))) 
    hIR_error = (fd.sqrt(fd.assemble(fd.dot(hIR - hI, hIR - hI) * fd.dx)))
    uTR_error = (fd.sqrt(fd.assemble(fd.dot(uTR - uI, uTR - uI) * fd.dx))) 
    hTR_error = (fd.sqrt(fd.assemble(fd.dot(hTR - hI, hTR - hI) * fd.dx)))
    uIE_error = (fd.sqrt(fd.assemble(fd.dot(uIE - uI, uIE - uI) * fd.dx))) 
    hIE_error = (fd.sqrt(fd.assemble(fd.dot(hIE - hI, hIE - hI) * fd.dx)))    
    uTE_error = (fd.sqrt(fd.assemble(fd.dot(uTE - uI, uTE - uI) * fd.dx))) 
    hTE_error = (fd.sqrt(fd.assemble(fd.dot(hTE - hI, hTE - hI) * fd.dx)))  

    if i == 1.0:
        u_ref = u_error
        h_ref = h_error
    #     uIR_ref = uIR_error
    #     hIR_ref = hIR_error
    #     uTR_ref = uTR_error
    #     hTR_ref = hTR_error
    #     uIE_ref = uIE_error
    #     hIE_ref = hIE_error
    #     uTE_ref = uTE_error
    #     hTE_ref = hTE_error        

    U = np.append(U, u_error)
    H = np.append(H, h_error)
    UIR = np.append(UIR, uIR_error)
    HIR = np.append(HIR, hIR_error)
    UTR = np.append(UTR, uTR_error)
    HTR = np.append(HTR, hTR_error)
    UIE = np.append(UIE, uIE_error)
    HIE = np.append(HIE, hIE_error)
    UTE = np.append(UTE, uTE_error)
    HTE = np.append(HTE, hTE_error)

    # print(f"L2 error at dt: {T} is u: {U[i]}, h: {H[i]}")

plt.figure()
plt.plot(dt, U, '-x', label='velocity TR-BDF2')
plt.plot(dt, H, '-x', label='height TR-BDF2')
plt.plot(dt, UIR, '-x', label='velocity IM(R)')
plt.plot(dt, HIR, '-x', label='height IM(R)')
plt.plot(dt, UTR, '-x', label='velocity TR-BDF2(R)')
plt.plot(dt, HTR, '-x', label='height TR-BDF2(R)')
plt.plot(dt, UIE, '-x', label='velocity IM(E)')
plt.plot(dt, HIE, '-x', label='height IM(E)')
plt.plot(dt, UTE, '-x', label='velocity TR-BDF2(E)')
plt.plot(dt, HTE, '-x', label='height TR-BDF2(E)')
plt.xscale('log', base=2)
plt.yscale('log')
plt.legend()
plt.xlabel('time step, dt')
plt.ylabel('Normalized difference')
plt.title("Comparison of the L2 difference between IM and other methods")
plt.show()
plt.savefig('convergence_'+str(dmax)+'.png')
        