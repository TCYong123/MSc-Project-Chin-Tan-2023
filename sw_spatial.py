import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import matplotlib.pyplot as plt
import argparse
import numpy as np

parser = argparse.ArgumentParser()
args = parser.parse_known_args()
args = args[0]

T = 60*60
dmax = 50.0

im_energy0 = np.loadtxt("im_energy0_"+str(T)+"_"+str(dmax)+".array")
im_energy1 = np.loadtxt("im_energy1_"+str(T)+"_"+str(dmax)+".array")
im_energy2 = np.loadtxt("im_energy2_"+str(T)+"_"+str(dmax)+".array")
im_energy3 = np.loadtxt("im_energy3_"+str(T)+"_"+str(dmax)+".array")
im_energy4 = np.loadtxt("im_energy4_"+str(T)+"_"+str(dmax)+".array")

n0 = np.arange(len(im_energy0))
n1 = np.arange(len(im_energy1))
n2 = np.arange(len(im_energy2))
n3 = np.arange(len(im_energy3))
n4 = np.arange(len(im_energy4))

plt.figure()
plt.semilogy(n0, im_energy1, label="IM_1")
plt.semilogy(n1, im_energy1, label="IM_2")
plt.semilogy(n2, im_energy2, label="IM_3")
plt.semilogy(n3, im_energy3, label="IM_4")
plt.semilogy(n4, im_energy4, label="IM_5")
plt.legend()
plt.xlabel("Time (hours)")
plt.ylabel("Total energy in system")
plt.title("Energy comparison for different spatial discretization")
plt.show()




        