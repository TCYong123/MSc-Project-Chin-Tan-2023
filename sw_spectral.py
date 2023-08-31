import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.fft as fft

dt = 3600.0
dmax = 50.0

im_arr = np.loadtxt("im_ht"+str(dt)+"_"+str(dmax)+".array")
tr_arr = np.loadtxt("tr_ht"+str(dt)+"_"+str(dmax)+".array")
imR_arr = np.loadtxt("imR_ht"+str(dt)+"_"+str(dmax)+".array")
trR_arr = np.loadtxt("trR_ht"+str(dt)+"_"+str(dmax)+".array")
imE_arr = np.loadtxt("imE_ht"+str(dt)+"_"+str(dmax)+".array")
trE_arr = np.loadtxt("trE_ht"+str(dt)+"_"+str(dmax)+".array")


tr = np.arange(len(tr_arr))

im_arr -= np.mean(im_arr)
tr_arr -= np.mean(tr_arr)
imR_arr -= np.mean(imR_arr)
trR_arr -= np.mean(trR_arr)
imE_arr -= np.mean(imE_arr)
trE_arr -= np.mean(trE_arr)

mask = 0.5 - 0.5 * np.cos(2*np.pi*tr/len(tr))

F1 = fft.fft(mask*im_arr)
F2 = fft.fft(mask*tr_arr)
F3 = fft.fft(mask*imR_arr)
F4 = fft.fft(mask*trR_arr)
F5 = fft.fft(mask*imE_arr)
F6 = fft.fft(mask*trE_arr)

plt.figure()
plt.semilogy(abs(F1)**2 ,label='IM')
plt.semilogy(abs(F2)**2 ,label='TR-BDF2')
plt.legend()
plt.xlabel('Frequency')
plt.ylabel('Amplitude')
plt.title('Spectral Analysis of depth (Standard)')
plt.xlim(0, 605)
plt.savefig("spectral_s_"+str(dt)+"_"+str(dmax)+".png")

plt.figure()
plt.semilogy(abs(F3)**2 ,label='IM(R)')
plt.semilogy(abs(F4)**2 ,label='TR-BDF2(R)')
plt.legend()
plt.xlabel('Frequency')
plt.ylabel('Amplitude')
plt.title('Spectral Analysis of depth (Rosenbrock)')
plt.xlim(0, 605)
plt.savefig("spectral_r_"+str(dt)+"_"+str(dmax)+".png")

plt.figure()
plt.semilogy(abs(F5)**2 ,label='IM(E)')
plt.semilogy(abs(F6)**2 ,label='TR-BDF2(E)')
plt.legend()
plt.xlabel('Frequency')
plt.ylabel('Amplitude')
plt.title('Spectral Analysis of depth (Eisenstat-Walker)')
plt.xlim(0, 605)
plt.savefig("spectral_ew_"+str(dt)+"_"+str(dmax)+".png")