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
parser.add_argument('--dmax', type=float, default=50.0)
parser.add_argument('--dt', type=int, default=0)
args = parser.parse_known_args()
args = args[0]
QI = np.array([])
QT = np.array([])
QIR = np.array([])
QTR = np.array([])
QIE = np.array([])
QTE = np.array([])
dmax = args.dmax

if args.w:
    dt = [2.0**(-args.dt)]
else:    
    dt = [2.0**(-n) for n in range(-1,4)]

for T in dt:
    if args.w:
        Mesh = sw_create_mesh.main(["--filename=vort"+str(dt), "--ref_level=5", "--dmax="+str(dmax)])
        sw_im.main(["--filename=vort"+str(dt), "--write=3", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_trbdf2_R.main(["--filename=vort"+str(dt), "--write=3", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_im_R.main(["--filename=vort"+str(dt), "--write=3", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_trbdf2.main(["--filename=vort"+str(dt), "--write=3", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_im_E.main(["--filename=vort"+str(dt), "--write=3", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)
        sw_trbdf2_E.main(["--filename=vort"+str(dt), "--write=3", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], Mesh)

    with fd.CheckpointFile("vorticity_dt"+str(T*60*60)+"_"+str(dmax)+".h5", 'r') as afile:
        mesh = afile.load_mesh("Mesh")      
        qI = afile.load_function(mesh, "q_outI")
        qT = afile.load_function(mesh, "q_outT")
        qIR = afile.load_function(mesh, "q_outIR")
        qTR = afile.load_function(mesh, "q_outTR")
        qIE = afile.load_function(mesh, "q_outIE")
        qTE = afile.load_function(mesh, "q_outTE")

    qI_error = (fd.sqrt(fd.assemble(fd.dot(qI-qI, qI-qI) * fd.dx))) 
    qT_error = (fd.sqrt(fd.assemble(fd.dot(qT-qI, qT-qI) * fd.dx))) 
    qIR_error = (fd.sqrt(fd.assemble(fd.dot(qIR-qI, qIR-qI) * fd.dx))) 
    qTR_error = (fd.sqrt(fd.assemble(fd.dot(qTR-qI, qTR-qI) * fd.dx))) 
    # qIE_error = (fd.sqrt(fd.assemble(fd.dot(qIE, qIE) * fd.dx))) 
    # qTE_error = (fd.sqrt(fd.assemble(fd.dot(qTE, qTE) * fd.dx))) 

    QI = np.append(QI, qI_error)
    QT = np.append(QT, qT_error)
    QIR = np.append(QIR, qIR_error)
    QTR = np.append(QTR, qTR_error)
    # QIE = np.append(QIE, qIE_error)
    # QTE = np.append(QTE, qTE_error)


plt.figure()
# plt.plot(dt, QI, '--x', label='IM')
plt.plot(dt, QT, '-.o', label='TR-BDF2')
plt.plot(dt, QIR, '--^', label='IM(R)')
plt.plot(dt, QTR, '-*', label='TR-BDF2(R)')
# plt.plot(dt, QIE, '-x', label='IM(E)')
# plt.plot(dt, QTE, '-x', label='TR-BDF2(E)')
plt.xscale('log', base=2)
# plt.yscale('log')
plt.legend()
plt.xlabel('time step size, dt')
plt.ylabel('L2 difference')
plt.title("Potential vorticity, dmax=50.0")
# plt.show()
plt.savefig('vorticity_new_'+str(dt)+'.png', bbox_inches="tight")
        
