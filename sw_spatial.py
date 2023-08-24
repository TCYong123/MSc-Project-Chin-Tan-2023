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
parser.add_argument('--ref', type=int, default=5)
args = parser.parse_known_args()
args = args[0]
dmax = args.dmax
dt = 2.0**(-args.dt)

T = dt*60*60
dmax = args.dmax
ref = args.ref

if args.w:
    mesh = sw_create_mesh.main(["--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(dt)])
    sw_im.main(["--energy=2", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_trbdf2.main(["--energy=2", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_im_R.main(["--energy=2", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_trbdf2_R.main(["--energy=2", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_im_E.main(["--energy=2", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_trbdf2_E.main(["--energy=2", "--ref_level="+str(ref), "--dmax="+str(dmax), "--dt="+str(dt)], mesh) 

# im_energy0 = np.loadtxt("im_energy0"+"_"+str(T)+"_"+str(dmax)+".array")
# im_energy1 = np.loadtxt("im_energy1"+"_"+str(T)+"_"+str(dmax)+".array")
# im_energy2 = np.loadtxt("im_energy2"+"_"+str(T)+"_"+str(dmax)+".array")
# im_energy3 = np.loadtxt("im_energy3"+"_"+str(T)+"_"+str(dmax)+".array")
# im_energy4 = np.loadtxt("im_energy4"+"_"+str(T)+"_"+str(dmax)+".array")

# imR_energy0 = np.loadtxt("imR_energy0"+"_"+str(T)+"_"+str(dmax)+".array")
# imR_energy1 = np.loadtxt("imR_energy1"+"_"+str(T)+"_"+str(dmax)+".array")
# imR_energy2 = np.loadtxt("imR_energy2"+"_"+str(T)+"_"+str(dmax)+".array")
# imR_energy3 = np.loadtxt("imR_energy3"+"_"+str(T)+"_"+str(dmax)+".array")
# imR_energy4 = np.loadtxt("imR_energy4"+"_"+str(T)+"_"+str(dmax)+".array")

# imE_energy0 = np.loadtxt("imE_energy0"+"_"+str(T)+"_"+str(dmax)+".array")
# imE_energy1 = np.loadtxt("imE_energy1"+"_"+str(T)+"_"+str(dmax)+".array")
# imE_energy2 = np.loadtxt("imE_energy2"+"_"+str(T)+"_"+str(dmax)+".array")
# imE_energy3 = np.loadtxt("imE_energy3"+"_"+str(T)+"_"+str(dmax)+".array")
# imE_energy4 = np.loadtxt("imE_energy4"+"_"+str(T)+"_"+str(dmax)+".array")

# tr_energy0 = np.loadtxt("tr_energy0"+"_"+str(T)+"_"+str(dmax)+".array")
# tr_energy1 = np.loadtxt("tr_energy1"+"_"+str(T)+"_"+str(dmax)+".array")
# tr_energy2 = np.loadtxt("tr_energy2"+"_"+str(T)+"_"+str(dmax)+".array")
# tr_energy3 = np.loadtxt("tr_energy3"+"_"+str(T)+"_"+str(dmax)+".array")
# tr_energy4 = np.loadtxt("tr_energy4"+"_"+str(T)+"_"+str(dmax)+".array")

# trR_energy0 = np.loadtxt("trR_energy0"+"_"+str(T)+"_"+str(dmax)+".array")
# trR_energy1 = np.loadtxt("trR_energy1"+"_"+str(T)+"_"+str(dmax)+".array")
# trR_energy2 = np.loadtxt("trR_energy2"+"_"+str(T)+"_"+str(dmax)+".array")
# trR_energy3 = np.loadtxt("trR_energy3"+"_"+str(T)+"_"+str(dmax)+".array")
# trR_energy4 = np.loadtxt("trR_energy4"+"_"+str(T)+"_"+str(dmax)+".array")

# trE_energy0 = np.loadtxt("trE_energy0"+"_"+str(T)+"_"+str(dmax)+".array")
# trE_energy1 = np.loadtxt("trE_energy1"+"_"+str(T)+"_"+str(dmax)+".array")
# trE_energy2 = np.loadtxt("trE_energy2"+"_"+str(T)+"_"+str(dmax)+".array")
# trE_energy3 = np.loadtxt("trE_energy3"+"_"+str(T)+"_"+str(dmax)+".array")
# trE_energy4 = np.loadtxt("trE_energy4"+"_"+str(T)+"_"+str(dmax)+".array")

# n0 = np.arange(len(im_energy0))
# n1 = np.arange(len(im_energy1))
# n2 = np.arange(len(im_energy2))
# n3 = np.arange(len(im_energy3))
# n4 = np.arange(len(im_energy4))

# plt.figure()
# plt.semilogy(n0, im_energy1, label="IM_1")
# plt.semilogy(n1, im_energy1, label="IM_2")
# plt.semilogy(n2, im_energy2, label="IM_3")
# plt.semilogy(n3, im_energy3, label="IM_4")
# plt.semilogy(n4, im_energy4, label="IM_5")
# plt.legend()
# plt.xlabel("Time (hours)")
# plt.ylabel("Total energy in system")
# plt.title("Energy comparison for different spatial discretization, dt=1.0")
# plt.savefig('spatial_im_1.0_50.png')

# n0 = np.arange(len(imR_energy0))
# n1 = np.arange(len(imR_energy1))
# n2 = np.arange(len(imR_energy2))
# n3 = np.arange(len(imR_energy3))
# n4 = np.arange(len(imR_energy4))

# plt.figure()
# plt.semilogy(n0, imR_energy1, label="IM(R)_1")
# plt.semilogy(n1, imR_energy1, label="IM(R)_2")
# plt.semilogy(n2, imR_energy2, label="IM(R)_3")
# plt.semilogy(n3, imR_energy3, label="IM(R)_4")
# plt.semilogy(n4, imR_energy4, label="IM(R)_5")
# plt.legend()
# plt.xlabel("Time (hours)")
# plt.ylabel("Total energy in system")
# plt.title("Energy comparison for different spatial discretization, dt=1.0")
# plt.savefig('spatial_imR_1.0_50.png')


# n0 = np.arange(len(imE_energy0))
# n1 = np.arange(len(imE_energy1))
# n2 = np.arange(len(imE_energy2))
# n3 = np.arange(len(imE_energy3))
# n4 = np.arange(len(imE_energy4))

# plt.figure()
# plt.semilogy(n0, imE_energy1, label="IM(E)_1")
# plt.semilogy(n1, imE_energy1, label="IM(E)_2")
# plt.semilogy(n2, imE_energy2, label="IM(E)_3")
# plt.semilogy(n3, imE_energy3, label="IM(E)_4")
# plt.semilogy(n4, imE_energy4, label="IM(E)_5")
# plt.legend()
# plt.xlabel("Time (hours)")
# plt.ylabel("Total energy in system")
# plt.title("Energy comparison for different spatial discretization, dt=1.0")
# plt.savefig('spatial_imE_1.0_50.png')

# n0 = np.arange(len(tr_energy0))
# n1 = np.arange(len(tr_energy1))
# n2 = np.arange(len(tr_energy2))
# n3 = np.arange(len(tr_energy3))
# n4 = np.arange(len(tr_energy4))

# plt.figure()
# plt.semilogy(n0, tr_energy1, label="TR-BDF2_1")
# plt.semilogy(n1, tr_energy1, label="TR-BDF2_2")
# plt.semilogy(n2, tr_energy2, label="TR-BDF2_3")
# plt.semilogy(n3, tr_energy3, label="TR-BDF2_4")
# plt.semilogy(n4, tr_energy4, label="TR-BDF2_5")
# plt.legend()
# plt.xlabel("Time (hours)")
# plt.ylabel("Total energy in system")
# plt.title("Energy comparison for different spatial discretization, dt=1.0")
# plt.savefig('spatial_tr_1.0_50.png')

# n0 = np.arange(len(trR_energy0))
# n1 = np.arange(len(trR_energy1))
# n2 = np.arange(len(trR_energy2))
# n3 = np.arange(len(trR_energy3))
# n4 = np.arange(len(trR_energy4))

# plt.figure()
# plt.semilogy(n0, trR_energy1, label="TR-BDF2(R)_1")
# plt.semilogy(n1, trR_energy1, label="TR-BDF2(R)_2")
# plt.semilogy(n2, trR_energy2, label="TR-BDF2(R)_3")
# plt.semilogy(n3, trR_energy3, label="TR-BDF2(R)_4")
# plt.semilogy(n4, trR_energy4, label="TR-BDF2(R)_5")
# plt.legend()
# plt.xlabel("Time (hours)")
# plt.ylabel("Total energy in system")
# plt.title("Energy comparison for different spatial discretization, dt=1.0")
# plt.savefig('spatial_trR_1.0_50.png')


# n0 = np.arange(len(trE_energy0))
# n1 = np.arange(len(trE_energy1))
# n2 = np.arange(len(trE_energy2))
# n3 = np.arange(len(trE_energy3))
# n4 = np.arange(len(trE_energy4))

# plt.figure()
# plt.semilogy(n0, trE_energy1, label="TR-BDF2(E)_1")
# plt.semilogy(n1, trE_energy1, label="TR-BDF2(E)_2")
# plt.semilogy(n2, trE_energy2, label="TR-BDF2(E)_3")
# plt.semilogy(n3, trE_energy3, label="TR-BDF2(E)_4")
# plt.semilogy(n4, trE_energy4, label="TR-BDF2(E)_5")
# plt.legend()
# plt.xlabel("Time (hours)")
# plt.ylabel("Total energy in system")
# plt.title("Energy comparison for different spatial discretization, dt=1.0")
# plt.savefig('spatial_trE_1.0_50.png')        