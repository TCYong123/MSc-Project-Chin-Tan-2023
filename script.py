import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
import firedrake as fd
import sw_im
import sw_trbdf2
import sw_create_mesh
import argparse

parser = argparse.ArgumentParser()
args = parser.parse_known_args()
args = args[0]

T = 1.0

sw_im.main(["--ref_level=2"])




        