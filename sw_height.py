import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import numpy as np
import matplotlib.pyplot as plt

dt = 3600.0

im_arr = np.loadtxt("im_ht"+str(dt)+"_50.0.array")
tr_arr = np.loadtxt("tr_ht"+str(dt)+"_50.0.array")
T = np.arange(0, len(im_arr))*dt/3600.0

plt.figure()
plt.plot(T, im_arr, label="im")
plt.legend()
plt.xlabel("Time in hours")
plt.ylabel("Height in meters")
plt.title("Height of a point against time (IM) "+str(dt))
plt.savefig("height_im"+str(dt)+"_50.png")

plt.figure()
plt.plot(T, tr_arr, label="trbdf2")
plt.legend()
plt.xlabel("Time in hours")
plt.ylabel("Height in meters")
plt.title("Height of a point against time (TR-BDF2) "+str(dt))
plt.savefig("height_tr"+str(dt)+"_50.png")

plt.figure()
plt.plot(T, im_arr, label="im")
plt.plot(T, tr_arr, label="trbdf2")
plt.legend()
plt.xlabel("Time in hours")
plt.ylabel("Height in meters")
plt.title("Height of a point against time "+str(dt))
plt.savefig("height"+str(dt)+"_50.png")