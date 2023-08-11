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
parser.add_argument('--ref', type=int, default=5)
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

if args.w:
    ref_level = [args.ref]
else:    
    ref_level = [1, 2, 3, 4, 5]

T = 1.0

for ref in ref_level:
    if args.w:
        Mesh = sw_create_mesh.main(["--ref_level="+str(ref), "--dmax="+str(dmax)])
        sw_im.main(["--write=21", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt=0.0009765625"], Mesh)
        sw_im.main(["--write=2", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_trbdf2_R.main(["--write=2", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_im_R.main(["--write=2", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_trbdf2.main(["--write=2", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_im_E.main(["--write=2", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_trbdf2_E.main(["--write=2", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(T)], Mesh)

    with fd.CheckpointFile("convergence_ds"+str(T*60*60)+"_"+str(ref)+"_"+str(dmax)+".h5", 'r') as afile:
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

    UI = np.append(U, uI_error)
    HI = np.append(H, hI_error)
    UT = np.append(U, uT_error)
    HT = np.append(H, hT_error)
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
plt.plot(ref_level, UI, '-x', label='velocity IM')
plt.plot(ref_level, HI, '-x', label='height IM')
plt.plot(ref_level, UT, '-x', label='velocity TR-BDF2')
plt.plot(ref_level, HT, '-x', label='height TR-BDF2')
plt.plot(ref_level, UIR, '-x', label='velocity IM(R)')
plt.plot(ref_level, HIR, '-x', label='height IM(R)')
plt.plot(ref_level, UTR, '-x', label='velocity TR-BDF2(R)')
plt.plot(ref_level, HTR, '-x', label='height TR-BDF2(R)')
plt.plot(ref_level, UIE, '-x', label='velocity IM(E)')
plt.plot(ref_level, HIE, '-x', label='height IM(E)')
plt.plot(ref_level, UTE, '-x', label='velocity TR-BDF2(E)')
plt.plot(ref_level, HTE, '-x', label='height TR-BDF2(E)')
plt.xscale('log', base=2)
plt.yscale('log')
plt.legend()
plt.xlabel('Spatial discretization level')
plt.ylabel('Relative L2 error')
plt.title("Comparison of the L2 difference of the methods")
plt.savefig('convergence_spatial_'+str(dmax)+'.png')
        