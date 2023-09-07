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
parser.add_argument('--dt', type=float, default=1)
args = parser.parse_known_args()
args = args[0]
dmax = args.dmax
dt = 1.0*(args.dt)

if args.w == True:
    mesh = sw_create_mesh.main(["--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)])
    sw_im.main(["--energy=1", "--iter=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_trbdf2.main(["--energy=1", "--iter=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_im_R.main(["--energy=1", "--iter=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_trbdf2_R.main(["--energy=1", "--iter=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_im_E.main(["--energy=1", "--iter=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)
    sw_trbdf2_E.main(["--energy=1", "--iter=1", "--ref_level=5", "--dmax="+str(dmax), "--dt="+str(dt)], mesh)

T = dt*60*60
im_energy = np.loadtxt("im_benergy"+str(T)+"_"+str(dmax)+".array") 
tr_energy = np.loadtxt("tr_benergy"+str(T)+"_"+str(dmax)+".array")
imR_energy = np.loadtxt("imR_benergy"+str(T)+"_"+str(dmax)+".array")
trR_energy= np.loadtxt("trR_benergy"+str(T)+"_"+str(dmax)+".array")
# imE_energy = np.loadtxt("imE_benergy"+str(T)+"_"+str(dmax)+".array")
# trE_energy= np.loadtxt("trE_benergy"+str(T)+"_"+str(dmax)+".array")

im_energy2 = (im_energy - im_energy[0]) /im_energy[0]
tr_energy2 = (tr_energy - tr_energy[0]) /tr_energy[0]
imR_energy2 = (imR_energy - imR_energy[0]) /imR_energy[0]
trR_energy2= (trR_energy - trR_energy[0]) /trR_energy[0]
# imE_energy2 = (imE_energy - imE_energy[0]) /imE_energy[0]
# trE_energy2= (trE_energy - trE_energy[0]) /trE_energy[0]

n = np.arange(len(im_energy)*dt, step=dt)

# plt.figure()
plt.plot(n, im_energy2,linewidth=2,linestyle='--', label="IM")
plt.plot(n, tr_energy2,linewidth=2, label="TR-BDF2")
plt.plot(n, imR_energy2,linewidth=2,linestyle='-.', label="IM(R)")
plt.plot(n, trR_energy2,linewidth=2, linestyle=':',label="TR-BDF2(R)")
# plt.plot(n, imE_energy2, label="IM(E)")
# plt.plot(n, trE_energy2, label="TR-BDF2(E)")
plt.legend()
plt.ylabel('Relative energy difference')
plt.xlabel("Time (hours)")
plt.title("Energy, dt="+str(dt)+", dmax="+str(dmax))
# plt.ylim((-3.5*10**-6,1.5*10**-6))
# plt.show()
plt.savefig('benergy_'+str(T)+'_'+str(dmax)+'.png', bbox_inches="tight")


# Equation for energy: "E = KE + PE = (h*|U|^2/2 + g*(h**2/2 + h*b)*dx"

im_itcount = np.loadtxt("im_itcount_"+str(T)+"_"+str(dmax)+".array")
tr_itcount = np.loadtxt("tr_itcount_"+str(T)+"_"+str(dmax)+".array")
imR_itcount = np.loadtxt("imR_itcount_"+str(T)+"_"+str(dmax)+".array")
trR_itcount= np.loadtxt("trR_itcount_"+str(T)+"_"+str(dmax)+".array")
imE_itcount = np.loadtxt("imE_itcount_"+str(T)+"_"+str(dmax)+".array")
trE_itcount= np.loadtxt("trE_itcount_"+str(T)+"_"+str(dmax)+".array")

im_stepcount = np.loadtxt("im_stepcount_"+str(T)+"_"+str(dmax)+".array")
tr_stepcount = np.loadtxt("tr_stepcount_"+str(T)+"_"+str(dmax)+".array")
imR_stepcount = np.loadtxt("imR_stepcount_"+str(T)+"_"+str(dmax)+".array")
trR_stepcount= np.loadtxt("trR_stepcount_"+str(T)+"_"+str(dmax)+".array")
imE_stepcount = np.loadtxt("imE_stepcount_"+str(T)+"_"+str(dmax)+".array")
trE_stepcount= np.loadtxt("trE_stepcount_"+str(T)+"_"+str(dmax)+".array")

plt.figure()
plt.plot(im_stepcount, im_itcount, '-',label='IM')
plt.plot(imR_stepcount, imR_itcount, '-.',label='IM(R)')
plt.plot(imE_stepcount, imE_itcount, '--',label='IM(E)')
plt.legend()
plt.title("Newton's solve, dt="+str(dt)+", dmax="+str(dmax))
plt.xlabel('Stepcount')
plt.ylabel('Total number of linear solves')          
plt.savefig('iteration_im_'+str(T)+'_'+str(dmax)+'.png')

plt.figure()
plt.plot(tr_stepcount, tr_itcount, 'b',linestyle='-',label='TR-BDF2')
plt.plot(trR_stepcount, trR_itcount, 'r',linestyle='-.',label='TR-BDF2(R)')
plt.plot(trE_stepcount, trE_itcount, 'g',linestyle='--',label='TR-BDF2(E)')
plt.legend()
plt.title("Newton's solve, dt="+str(dt)+", dmax="+str(dmax))
plt.xlabel('Stepcount')
plt.ylabel('Total number of linear solves')          
# plt.show()
plt.savefig('iteration_tr_'+str(T)+'_'+str(dmax)+'.png')





