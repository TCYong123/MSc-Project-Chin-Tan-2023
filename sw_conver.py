import firedrake as fd
import sw_test
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--w', type=bool, default=False) #Toggle writing of files
args = parser.parse_known_args()
args = args[0]
U = []
H = []

for i, T in enumerate([2.0**(-n) for n in range(-2, 6)]):
    if args.w:
        sw_test.main(["--ref_level=3", "--dmax=1", "--dt="+str(T)])
    
    with fd.CheckpointFile("test_dt"+str(60*60*T)+".h5", 'r') as afile:
        meshI = afile.load_mesh("mesh")
        uI = afile.load_function(meshI, "uI")
        hI = afile.load_function(meshI, "hI")
        uT = afile.load_function(meshI, "uT")
        hT = afile.load_function(meshI, "hT")

    u_error = (fd.sqrt(fd.assemble(fd.dot(uT - uI, uT - uI) * fd.dx)))
    h_error = (fd.sqrt(fd.assemble(fd.dot(hT - hI, hT - hI) * fd.dx)))

    U.append(u_error)
    H.append(h_error)
    print(f"L2 error at dt: {T} is u: {U[i]}, h: {H[i]}")