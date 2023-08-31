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
args = parser.parse_known_args()
args = args[0]
dmax = args.dmax
dt = 2.0**(-args.dt)

if args.w == True:
    mesh = sw_create_mesh.main(["--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)])
    sw_im.main(["--array=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_trbdf2.main(["--array=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_im_R.main(["--array=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_trbdf2_R.main(["--array=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_im_E.main(["--array=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_trbdf2_E.main(["--array=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)

dt = 60*60*dt
im_arr = np.loadtxt("im_ht"+str(dt)+"_"+str(dmax)+".array")
tr_arr = np.loadtxt("tr_ht"+str(dt)+"_"+str(dmax)+".array")
imR_arr = np.loadtxt("imR_ht"+str(dt)+"_"+str(dmax)+".array")
trR_arr = np.loadtxt("trR_ht"+str(dt)+"_"+str(dmax)+".array")
imE_arr = np.loadtxt("imE_ht"+str(dt)+"_"+str(dmax)+".array")
trE_arr = np.loadtxt("trE_ht"+str(dt)+"_"+str(dmax)+".array")
T = np.arange(0, len(im_arr))*dt/3600.0

# plt.figure()
# plt.plot(T, im_arr, label="im")
# plt.legend()
# plt.xlabel("Time in hours")
# plt.ylabel("Height in metres")
# plt.title("Height of a point against time (IM) "+str(dt))
# plt.savefig("height_im"+str(dt)+"_50.png")

# plt.figure()
# plt.plot(T, tr_arr, label="trbdf2")
# plt.legend()
# plt.xlabel("Time in hours")
# plt.ylabel("Height in metres")
# plt.title("Height of a point against time (TR-BDF2) "+str(dt))
# plt.savefig("height_tr"+str(dt)+"_50.png")

plt.figure()
plt.plot(T, im_arr,label="IM")
plt.plot(T, tr_arr, label="TR-BDF2")
plt.legend()
plt.xlabel("Time (hours)")
plt.ylabel("Depth (metres)")
plt.title('Variation of height (Standard)')
plt.savefig("height_s"+str(dt)+"_"+str(dmax)+".png")

plt.figure()
plt.plot(T, imR_arr, label="IM(R)")
plt.plot(T, trR_arr, label="TR-BDF2(R)")
plt.legend()
plt.xlabel("Time (hours)")
plt.ylabel("Depth (metres)")
plt.title('Variation of height (Rosenbrock)')
plt.savefig("height_r"+str(dt)+"_"+str(dmax)+".png")

plt.figure()
plt.plot(T, imE_arr, label="IM(E)")
plt.plot(T, trE_arr, label="TR-BDF2(E)")
plt.legend()
plt.xlabel("Time (hours)")
plt.ylabel("Depth (metres)")
plt.title('Variation of height (Eisenstat-Walker)')
plt.savefig("height_ew"+str(dt)+"_"+str(dmax)+".png")