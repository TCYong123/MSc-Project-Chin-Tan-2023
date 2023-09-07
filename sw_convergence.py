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
parser.add_argument('--dt', type=int, default=1)
args = parser.parse_known_args()
args = args[0]
UI = np.array([])
HI = np.array([])
UT = np.array([])
HT = np.array([])
UIR = np.array([])
HIR = np.array([])
UTR = np.array([])
HTR = np.array([])
UIE = np.array([])
HIE = np.array([])
UTE = np.array([])
HTE = np.array([])
dmax = args.dmax

if args.w:
    dt = [2.0**(-args.dt)]
else:    
    dt = [2.0**(-n) for n in range(-2,9)]

for T in dt:
    if args.w:
        Mesh = sw_create_mesh.main(["--ref_level=2", "--dmax="+str(dmax)])
        sw_im.main(["--write=11", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_im.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_trbdf2_R.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_im_R.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_trbdf2.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_im_E.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_trbdf2_E.main(["--write=1", "--ref_level=2", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)

    with fd.CheckpointFile("convergence_dt"+str(T*60*60)+"_"+str(dmax)+".h5", 'r') as afile:
        mesh = afile.load_mesh("Mesh")
        u_ref = afile.load_function(mesh, "u_outref")
        h_ref = afile.load_function(mesh, "h_outref")        
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


    uI_error = (fd.sqrt(fd.assemble(fd.dot(uI - u_ref, uI - u_ref) * fd.dx))) 
    hI_error = (fd.sqrt(fd.assemble(fd.dot(hI - h_ref, hI - h_ref) * fd.dx)))
    uT_error = (fd.sqrt(fd.assemble(fd.dot(uT - u_ref, uT - u_ref) * fd.dx))) 
    hT_error = (fd.sqrt(fd.assemble(fd.dot(hT - h_ref, hT - h_ref) * fd.dx)))
    uIR_error = (fd.sqrt(fd.assemble(fd.dot(uIR - u_ref, uIR - u_ref) * fd.dx))) 
    hIR_error = (fd.sqrt(fd.assemble(fd.dot(hIR - h_ref, hIR - h_ref) * fd.dx)))
    uTR_error = (fd.sqrt(fd.assemble(fd.dot(uTR - u_ref, uTR - u_ref) * fd.dx))) 
    hTR_error = (fd.sqrt(fd.assemble(fd.dot(hTR - h_ref, hTR - h_ref) * fd.dx)))
    uIE_error = (fd.sqrt(fd.assemble(fd.dot(uIE - u_ref, uIE - u_ref) * fd.dx))) 
    hIE_error = (fd.sqrt(fd.assemble(fd.dot(hIE - h_ref, hIE - h_ref) * fd.dx)))    
    uTE_error = (fd.sqrt(fd.assemble(fd.dot(uTE - u_ref, uTE - u_ref) * fd.dx))) 
    hTE_error = (fd.sqrt(fd.assemble(fd.dot(hTE - h_ref, hTE - h_ref) * fd.dx)))       

    UI = np.append(UI, uI_error)
    HI = np.append(HI, hI_error)
    UT = np.append(UT, uT_error)
    HT = np.append(HT, hT_error)
    UIR = np.append(UIR, uIR_error)
    HIR = np.append(HIR, hIR_error)
    UTR = np.append(UTR, uTR_error)
    HTR = np.append(HTR, hTR_error)
    UIE = np.append(UIE, uIE_error)
    HIE = np.append(HIE, hIE_error)
    UTE = np.append(UTE, uTE_error)
    HTE = np.append(HTE, hTE_error)

    # print(f"L2 error at dt: {T} is u: {U[i]}, h: {H[i]}")
np.savetxt('UI.array',UI)
np.savetxt('UIR.array',UIR)
np.savetxt('UIE.array',UIE)
np.savetxt('UT.array',UT)
np.savetxt('UTR.array',UTR)
np.savetxt('UTE.array',UTE)

np.savetxt('HI.array',HI)
np.savetxt('HIR.array',HIR)
np.savetxt('HIE.array',HIE)
np.savetxt('HT.array',HT)
np.savetxt('HTR.array',HTR)
np.savetxt('HTE.array',HTE)

plt.figure()
plt.plot(dt, UI, '--X', label='velocity IM')
plt.plot(dt, HI, '-.^', label='height IM')
plt.plot(dt, UT, ':*', label='velocity TR-BDF2')
plt.plot(dt, HT, '-p', label='height TR-BDF2')
plt.xscale('log', base=2)
plt.yscale('log')
plt.legend()
plt.xlabel('time step size, dt')
plt.ylabel('Relative L2 error')
plt.title("Depth and velocity (Standard)")

plt.savefig('convergence_new_s'+str(dmax)+'.png')

plt.figure()
plt.plot(dt, UIR, '--X', label='velocity IM(R)')
plt.plot(dt, HIR, '-.^', label='height IM(R)')
plt.plot(dt, UTR, ':*', label='velocity TR-BDF2(R)')
plt.plot(dt, HTR, '-p', label='height TR-BDF2(R)')
plt.xscale('log', base=2)
plt.yscale('log')
plt.legend()
plt.xlabel('time step size, dt')
plt.ylabel('Relative L2 error')
plt.title("Depth and velocity (Rosenbrock)")

plt.savefig('convergence_new_r'+str(dmax)+'.png')

# plt.figure()
# plt.plot(dt, UIE, '--X', label='velocity IM(E)')
# plt.plot(dt, HIE, '-.^', label='height IM(E)')
# plt.plot(dt, UTE, ':*', label='velocity TR-BDF2(E)')
# plt.plot(dt, HTE, '-p', label='height TR-BDF2(E)')
# plt.xscale('log', base=2)
# plt.yscale('log')
# plt.legend()
# plt.xlabel('time step size, dt')
# plt.ylabel('Relative L2 error')
# plt.title("Comparison of L2 difference (Eisenstat-Walker)")

# plt.savefig('convergence_new_ew'+str(dmax)+'.png')
        