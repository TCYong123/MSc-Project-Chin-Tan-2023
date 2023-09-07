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

im_energy0 = np.loadtxt("im_benergy0"+str(T)+"_"+str(dmax)+".array")
im_energy1 = np.loadtxt("im_benergy1"+str(T)+"_"+str(dmax)+".array")
im_energy2 = np.loadtxt("im_benergy2"+str(T)+"_"+str(dmax)+".array")
im_energy3 = np.loadtxt("im_benergy3"+str(T)+"_"+str(dmax)+".array")
im_energy4 = np.loadtxt("im_benergy4"+str(T)+"_"+str(dmax)+".array")

imR_energy0 = np.loadtxt("imR_benergy0"+"_"+str(T)+"_"+str(dmax)+".array")
imR_energy1 = np.loadtxt("imR_benergy1"+"_"+str(T)+"_"+str(dmax)+".array")
imR_energy2 = np.loadtxt("imR_benergy2"+"_"+str(T)+"_"+str(dmax)+".array")
imR_energy3 = np.loadtxt("imR_benergy3"+"_"+str(T)+"_"+str(dmax)+".array")
imR_energy4 = np.loadtxt("imR_benergy4_"+str(T)+"_"+str(dmax)+".array")

# imE_energy0 = np.loadtxt("imE_benergy0"+str(T)+"_"+str(dmax)+".array")
# imE_energy1 = np.loadtxt("imE_benergy1"+str(T)+"_"+str(dmax)+".array")
# imE_energy2 = np.loadtxt("imE_benergy2"+str(T)+"_"+str(dmax)+".array")
# imE_energy3 = np.loadtxt("imE_benergy3"+str(T)+"_"+str(dmax)+".array")
# imE_energy4 = np.loadtxt("imE_benergy4"+str(T)+"_"+str(dmax)+".array")

tr_energy0 = np.loadtxt("tr_benergy0"+"_"+str(T)+"_"+str(dmax)+".array")
tr_energy1 = np.loadtxt("tr_benergy1"+"_"+str(T)+"_"+str(dmax)+".array")
tr_energy2 = np.loadtxt("tr_benergy2"+"_"+str(T)+"_"+str(dmax)+".array")
tr_energy3 = np.loadtxt("tr_benergy3"+"_"+str(T)+"_"+str(dmax)+".array")
tr_energy4 = np.loadtxt("tr_benergy4_"+str(T)+"_"+str(dmax)+".array")

trR_energy0 = np.loadtxt("trR_benergy0"+"_"+str(T)+"_"+str(dmax)+".array")
trR_energy1 = np.loadtxt("trR_benergy1"+"_"+str(T)+"_"+str(dmax)+".array")
trR_energy2 = np.loadtxt("trR_benergy2"+"_"+str(T)+"_"+str(dmax)+".array")
trR_energy3 = np.loadtxt("trR_benergy3"+"_"+str(T)+"_"+str(dmax)+".array")
trR_energy4 = np.loadtxt("trR_benergy4_"+str(T)+"_"+str(dmax)+".array")

# trE_energy0 = np.loadtxt("trE_benergy0"+"_"+str(T)+"_"+str(dmax)+".array")
# trE_energy1 = np.loadtxt("trE_benergy1"+"_"+str(T)+"_"+str(dmax)+".array")
# trE_energy2 = np.loadtxt("trE_benergy2"+"_"+str(T)+"_"+str(dmax)+".array")
# trE_energy3 = np.loadtxt("trE_benergy3"+"_"+str(T)+"_"+str(dmax)+".array")
# trE_energy4 = np.loadtxt("trE_benergy"+str(T)+"_"+str(dmax)+".array")


im_energy0 = np.abs(im_energy0 - im_energy0[0]) /im_energy0[0]
im_energy1 = np.abs(im_energy1 - im_energy1[0]) /im_energy1[0]
im_energy2 = np.abs(im_energy2 - im_energy2[0]) /im_energy2[0]
im_energy3 = np.abs(im_energy3 - im_energy3[0]) /im_energy3[0]
im_energy4 = np.abs(im_energy4 - im_energy4[0]) /im_energy4[0]

imR_energy0 = np.abs(imR_energy0 - imR_energy0[0]) /imR_energy0[0]
imR_energy1 = np.abs(imR_energy1 - imR_energy1[0]) /imR_energy1[0]
imR_energy2 = np.abs(imR_energy2 - imR_energy2[0]) /imR_energy2[0]
imR_energy3 = np.abs(imR_energy3 - imR_energy3[0]) /imR_energy3[0]
imR_energy4 = np.abs(imR_energy4 - imR_energy4[0]) /imR_energy4[0]

# imE_energy0 = imE_energy0/imE_energy0[0]
# imE_energy1 = imE_energy1/imE_energy1[0]
# imE_energy2 = imE_energy2/imE_energy2[0]
# imE_energy3 = imE_energy3/imE_energy3[0]
# imE_energy4 = imE_energy4/imE_energy4[0]

tr_energy0 = np.abs(tr_energy0 - tr_energy0[0]) /tr_energy0[0]
tr_energy1 = np.abs(tr_energy1 - tr_energy1[0]) /tr_energy1[0]
tr_energy2 = np.abs(tr_energy2 - tr_energy2[0]) /tr_energy2[0]
tr_energy3 = np.abs(tr_energy3 - tr_energy3[0]) /tr_energy3[0]
tr_energy4 = np.abs(tr_energy4 - tr_energy4[0]) /tr_energy4[0]

trR_energy0 = np.abs(trR_energy0 - trR_energy0[0]) /trR_energy0[0]
trR_energy1 = np.abs(trR_energy1 - trR_energy1[0]) /trR_energy1[0]
trR_energy2 = np.abs(trR_energy2 - trR_energy2[0]) /trR_energy2[0]
trR_energy3 = np.abs(trR_energy3 - trR_energy3[0]) /trR_energy3[0]
trR_energy4 = np.abs(trR_energy4 - trR_energy4[0]) /trR_energy4[0]

# trE_energy0 = trE_energy0/trE_energy0[0]
# trE_energy1 = trE_energy1/trE_energy1[0]
# trE_energy2 = trE_energy2/trE_energy2[0]
# trE_energy3 = trE_energy3/trE_energy3[0]
# trE_energy4 = trE_energy4/trE_energy4[0]

n0 = np.arange(len(im_energy0))
n1 = np.arange(len(im_energy1))
n2 = np.arange(len(im_energy2))
n3 = np.arange(len(im_energy3))
n4 = np.arange(len(im_energy4))
plt.figure()
plt.semilogy(n0, im_energy0, label="IM, ref=1")
plt.semilogy(n1, im_energy1, label="IM, ref=2")
plt.semilogy(n2, im_energy2, label="IM, ref=3")
plt.semilogy(n3, im_energy3, label="IM, ref=4")
plt.semilogy(n4, im_energy4, label="IM, ref=5")
plt.legend( loc='lower right')
plt.xlabel("Time (hours)")
plt.ylabel("Relative error for energy")
plt.title("Energy, dt=1.0, dmax=50.0")
# plt.ylim(0.99999)
# plt.show()
plt.savefig('aspatial_im_1.0_50.png', bbox_inches="tight")

n0 = np.arange(len(imR_energy0))
n1 = np.arange(len(imR_energy1))
n2 = np.arange(len(imR_energy2))
n3 = np.arange(len(imR_energy3))
n4 = np.arange(len(imR_energy4))
plt.figure()
plt.semilogy(n0, imR_energy0, label="IM(R), ref=1")
plt.semilogy(n1, imR_energy1, label="IM(R), ref=2")
plt.semilogy(n2, imR_energy2, label="IM(R), ref=3")
plt.semilogy(n3, imR_energy3, label="IM(R), ref=4")
plt.semilogy(n4, imR_energy4, label="IM(R), ref=5")
plt.legend( loc='lower right')
plt.xlabel("Time (hours)")
plt.ylabel("Relative error for energy")
plt.title("Energy, dt=1.0, dmax=50.0")
# plt.ylim(0.99999)
# plt.show()
plt.savefig('aspatial_imR_1.0_50.png', bbox_inches="tight")


# n0 = np.arange(len(imE_energy0))
# n1 = np.arange(len(imE_energy1))
# n2 = np.arange(len(imE_energy2))
# n3 = np.arange(len(imE_energy3))
# n4 = np.arange(len(imE_energy4))
# plt.figure()
# plt.semilogy(n0, imE_energy0, label="IM(E), ref=1")
# plt.semilogy(n1, imE_energy1, label="IM(E), ref=2")
# plt.semilogy(n2, imE_energy2, label="IM(E), ref=3")
# plt.semilogy(n3, imE_energy3, label="IM(E), ref=4")
# plt.semilogy(n4, imE_energy4, label="IM(E), ref=5")
# plt.legend( loc='lower right')
# plt.xlabel("Time (hours)")
# plt.ylabel("Normalized energy")
# plt.title("Energy comparison, dt=1.0, dmax=50.0")
# plt.ylim(0.99999)
# plt.show()
# plt.savefig('aspatial_imE_1.0_50.png', bbox_inches="tight")

n0 = np.arange(len(tr_energy0))
n1 = np.arange(len(tr_energy1))
n2 = np.arange(len(tr_energy2))
n3 = np.arange(len(tr_energy3))
n4 = np.arange(len(tr_energy4))
plt.figure()
plt.semilogy(n0, tr_energy0, label="TR-BDF2, ref=1")
plt.semilogy(n1, tr_energy1, label="TR-BDF2, ref=2")
plt.semilogy(n2, tr_energy2, label="TR-BDF2, ref=3")
plt.semilogy(n3, tr_energy3, label="TR-BDF2, ref=4")
plt.semilogy(n4, tr_energy4, label="TR-BDF2, ref=5")
plt.legend( loc='lower right')
plt.xlabel("Time (hours)")
plt.ylabel("Relative error for energy")
plt.title("Energy, dt=1.0, dmax=50.0")
# plt.ylim(0.99999)
# plt.show()
plt.savefig('aspatial_tr_1.0_50.png', bbox_inches="tight")

n0 = np.arange(len(trR_energy0))
n1 = np.arange(len(trR_energy1))
n2 = np.arange(len(trR_energy2))
n3 = np.arange(len(trR_energy3))
n4 = np.arange(len(trR_energy4))
plt.figure()
plt.semilogy(n0, trR_energy0, label="TR-BDF2(R), ref=1")
plt.semilogy(n1, trR_energy1, label="TR-BDF2(R), ref=2")
plt.semilogy(n2, trR_energy2, label="TR-BDF2(R), ref=3")
plt.semilogy(n3, trR_energy3, label="TR-BDF2(R), ref=4")
plt.semilogy(n4, trR_energy4, label="TR-BDF2(R), ref=5")
plt.legend( loc='lower right')
plt.xlabel("Time (hours)")
plt.ylabel("Relative error for energy")
plt.title("Energy, dt=1.0, dmax=50.0")
# plt.ylim(0.99999)
# plt.show()
plt.savefig('aspatial_trR_1.0_50.png', bbox_inches="tight")


# n0 = np.arange(len(trE_energy0))
# n1 = np.arange(len(trE_energy1))
# n2 = np.arange(len(trE_energy2))
# n3 = np.arange(len(trE_energy3))
# n4 = np.arange(len(trE_energy4))
# plt.figure()
# plt.semilogy(n0, trE_energy0, label="TR-BDF2(E), ref=1")
# plt.semilogy(n1, trE_energy1, label="TR-BDF2(E), ref=2")
# plt.semilogy(n2, trE_energy2, label="TR-BDF2(E), ref=3")
# plt.semilogy(n3, trE_energy3, label="TR-BDF2(E), ref=4")
# plt.semilogy(n4, trE_energy4, label="TR-BDF2(E), ref=5")
# plt.legend( loc='lower right')
# plt.xlabel("Time (hours)")
# plt.ylabel("Normalized energy")
# plt.ylim(0.99999)
# plt.title("Energy comparison, dt=1.0, dmax=50.0")
# plt.savefig('aspatial_trE_1.0_50.png', bbox_inches="tight")        