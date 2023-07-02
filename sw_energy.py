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
parser.add_argument('--dmax', type=float, default=15.0)
parser.add_argument('--dt', type=float, default=1.0)
args = parser.parse_known_args()
args = args[0]
dmax = args.dmax
T = args.dt

if args.w == True:
    mesh = sw_create_mesh.main(["--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)])
    sw_im.main(["--array=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
    sw_trbdf2.main(["--array=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
    sw_im_R.main(["--array=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], mesh)
    sw_trbdf2_R.main(["--array=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(T)], mesh)

im_height = np.loadtxt("im_ht"+str(T)+"_"+str(dmax)+".array")
imR_height = np.loadtxt("imR_ht"+str(T)+"_"+str(dmax)+".array")
tr_height = np.loadtxt("tr_ht"+str(T)+"_"+str(dmax)+".array")
trR_height = np.loadtxt("trR_ht"+str(T)+"_"+str(dmax)+".array")

im_b = np.loadtxt("im_b"+str(T)+"_"+str(dmax)+".array")
imR_b = np.loadtxt("imR_b"+str(T)+"_"+str(dmax)+".array")
tr_b = np.loadtxt("tr_b"+str(T)+"_"+str(dmax)+".array")
trR_b = np.loadtxt("trR_b"+str(T)+"_"+str(dmax)+".array")

im_velocity = np.loadtxt("im_vt"+str(T)+"_"+str(dmax)+".array")
imR_velocity = np.loadtxt("imR_vt"+str(T)+"_"+str(dmax)+".array")
tr_velocity = np.loadtxt("tr_vt"+str(T)+"_"+str(dmax)+".array")
trR_velocity = np.loadtxt("trR_vt"+str(T)+"_"+str(dmax)+".array")

g = fd.Constant(9.8)
dx = fd.dx

"E = KE + PE = (h*|U|^2 + g*(h+b))/2*dx"