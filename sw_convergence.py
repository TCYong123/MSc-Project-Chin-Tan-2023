import firedrake as fd
import sw_im
import sw_trbdf2
import sys

U = []
H = []

for T in [0.5, 1.0, 2.0, 4.0]:
    sw_im.main(["--ref_level=2", "--dmax=1", "--dt="+str(T)])
    sw_trbdf2.main(["--ref_level=2", "--dmax=1", "--dt="+str(T)])

    with fd.CheckpointFile("im_dt"+str(60*60*T)+".h5", 'r') as afile:
        meshI = afile.load_mesh("meshI")
        uI = afile.load_function(meshI, "u_out")
        hI = afile.load_function(meshI, "h_out")

    with fd.CheckpointFile("trbdf2_dt"+str(60*60*T)+".h5", 'r') as afile:
        meshT = afile.load_mesh("meshT")
        uT = afile.load_function(meshT, "u_out")
        hT = afile.load_function(meshT, "h_out")

    u_error = (fd.sqrt(fd.assemble(fd.dot(uT - uI, uT - uI) * fd.dx)))
    h_error = (fd.sqrt(fd.assemble(fd.dot(hT - hI, hT - hI) * fd.dx)))

    U.append(u_error)
    H.append(h_error)

for i, T in enumerate([0.5, 1, 5, 10]):
    print(f"L2 error at time {T} is u: {U[i]}, h: {H[i]}")