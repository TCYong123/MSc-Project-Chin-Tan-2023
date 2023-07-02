import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import numpy as np
import matplotlib.pyplot as plt

dt = 3600

im_arr = np.loadtxt("im_vor"+str(dt)+"_15.array")
tr_arr = np.loadtxt("tr_vor"+str(dt)+"_15.array")
T = np.arange(0, len(im_arr))*dt/3600.0


plt.figure()
plt.plot(T, im_arr, label="im")
plt.plot(T, tr_arr, label="trbdf2")
plt.legend()
plt.xlabel("Time in hours")
plt.ylabel("Vorticity")
plt.title("Vorticity of a point against time "+str(dt))
plt.savefig("vort"+str(dt)+"_15.png")