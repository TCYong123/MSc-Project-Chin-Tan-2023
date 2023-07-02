import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import numpy as np
import matplotlib.pyplot as plt

im_arr = np.loadtxt("im_ht3600.0_50.0.array")
tr_arr = np.loadtxt("trbdf2_ht3600.0_50.0.array")

im_arr = im_arr[600:]
tr_arr = tr_arr[600:]


T = np.arange(0, len(im_arr))

fft_im = np.fft.rfft(im_arr)
fft_tr = np.fft.rfft(tr_arr)    
N = len(fft_im)
fft_im = np.fft.fftshift(fft_im)/N
fft_tr = np.fft.fftshift(fft_tr)/N
n = np.arange(-N/2, N/2)
freq = n/(2*(np.pi))

plt.figure()
plt.semilogy(freq,np.abs(fft_im),label="IM")
plt.xlabel("Frequency")
plt.ylabel("Amplitude in log-scale")
plt.legend()
plt.title("Fourier analysis of height against time (IM)")
plt.xlim(-0.2,15)
plt.savefig("fourier50_im.png")

plt.figure()
plt.semilogy(freq,np.abs(fft_tr),label="TRBDF2")
plt.xlabel("Frequency")
plt.ylabel("Amplitude in log-scale")
plt.legend()
plt.title("Fourier analysis of height against time (TRBDF2)")
plt.xlim(-0.2,15)
plt.savefig("fourier50_tr.png")

plt.figure()
plt.semilogy(freq,np.abs(fft_im),label="IM")
plt.semilogy(freq,np.abs(fft_tr),label="TRBDF2")
plt.xlabel("Frequency")
plt.ylabel("Amplitude in log-scale")
plt.legend()
plt.title("Fourier analysis of height against time")
plt.xlim(-0.2,15)
plt.savefig("fourier50.png")