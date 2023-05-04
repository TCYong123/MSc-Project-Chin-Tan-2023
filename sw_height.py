import numpy as np
import matplotlib.pyplot as plt

im_arr = np.loadtxt("im_ht.array")
tr_arr = np.loadtxt("trbdf2_ht.array")
T = np.arange(0, len(im_arr))

plt.figure()
plt.plot(T, im_arr, label="im")
plt.legend()
plt.xlabel("Time in hours")
plt.ylabel("Height in meters")
plt.title("Height of a point against time (IM)")
plt.savefig("height_im.png")

plt.figure()
plt.plot(T, tr_arr, label="trbdf2")
plt.legend()
plt.xlabel("Time in hours")
plt.ylabel("Height in meters")
plt.title("Height of a point against time (TR-BDF2)")
plt.savefig("height_tr.png")

plt.figure()
plt.plot(T, im_arr, label="im")
plt.plot(T, tr_arr, label="trbdf2")
plt.legend()
plt.xlabel("Time in hours")
plt.ylabel("Height in meters")
plt.title("Height of a point against time")
plt.savefig("height.png")