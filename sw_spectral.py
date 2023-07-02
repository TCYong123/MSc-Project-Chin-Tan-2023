import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.fft as fft

dt = 14400.0

tr_arr = np.loadtxt("tr_ht"+str(dt)+"_50.0.array")
im_arr = np.loadtxt("im_ht"+str(dt)+"_50.0.array")

tr = np.arange(len(tr_arr))

tr_arr -= np.mean(tr_arr)
im_arr -= np.mean(im_arr)

mask = 0.5 - 0.5 * np.cos(2*np.pi*tr/len(tr))

F1 = fft.fft(mask*tr_arr)
F2 = fft.fft(mask*im_arr)


plt.semilogy(abs(F1)**2 ,label='tr')
plt.semilogy(abs(F2)**2 ,label='im')

plt.title(str(dt))
plt.savefig("spectral_"+str(dt)+".png")