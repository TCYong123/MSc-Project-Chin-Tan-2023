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
QT = np.array([])
QIR = np.array([])
QTR = np.array([])
QIE = np.array([])
QTE = np.array([])
dmax = args.dmax

if args.w:
    dt = [2.0**(-args.dt)]
else:    
    dt = [2.0**(-n) for n in range(-1,7)]

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

    qT_error = (fd.sqrt(fd.assemble(fd.dot(qT - qI, qT - qI) * fd.dx))) 
    qIR_error = (fd.sqrt(fd.assemble(fd.dot(qIR - qI, qIR - qI) * fd.dx))) 
    qTR_error = (fd.sqrt(fd.assemble(fd.dot(qTR - qI, qTR - qI) * fd.dx))) 
    qIE_error = (fd.sqrt(fd.assemble(fd.dot(qIE - qI, qIE - qI) * fd.dx))) 
    qTE_error = (fd.sqrt(fd.assemble(fd.dot(qTE - qI, qTE - qI) * fd.dx))) 

    QT = np.append(QT, qT_error)
    QIR = np.append(QIR, qIR_error)
    QTR = np.append(QTR, qTR_error)
    QIE = np.append(QIE, qIE_error)
    QTE = np.append(QTE, qTE_error)


plt.figure()
plt.plot(dt, QT, '-x', label='vorticity TR-BDF2')
plt.plot(dt, QIR, '-x', label='vorticity IM(R)')
plt.plot(dt, QTR, '-x', label='vorticity TR-BDF2(R)')
plt.plot(dt, QIE, '-x', label='vorticity IM(E)')
plt.plot(dt, QTE, '-x', label='vorticity TR-BDF2(E)')
plt.xscale('log', base=2)
plt.yscale('log')
plt.legend()
plt.xlabel('time step, dt')
plt.ylabel('Relative L2 error')
plt.title("Comparison of the L2 difference of the methods")
plt.show()
plt.savefig('vorticity_new_'+str(dt)+'.png')
        
