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


plt.semilogy(abs(F1)**2 ,label='IM')
plt.semilogy(abs(F2)**2 ,label='TR-BDF2')
plt.semilogy(abs(F3)**2 ,label='IM(R)')
plt.semilogy(abs(F4)**2 ,label='TR-BDF2(R)')
plt.semilogy(abs(F5)**2 ,label='IM(E)')
plt.semilogy(abs(F6)**2 ,label='TR-BDF2(E)')
plt.legend()
plt.xlabel('Frequency')
plt.ylabel('Amplitude')
plt.title('Spectral Analysis of depth, dt='+str(dt)+', dmax='+str(dmax))
plt.xlim(0, 605)
plt.savefig("spectral_"+str(dt)+"_"+str(dmax)+".png")