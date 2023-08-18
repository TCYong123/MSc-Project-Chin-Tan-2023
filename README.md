# MSc-Project-Chin-Tan-2023 
Last modified date: 18/08/2023

Main Files: \
sw_im.py: Crank-Nicholson and Implicit Midpoint methods \
sw_trbdf2.py: TR-BDF2 (Trapezium Rule - Backward Differentiation Order 2) method \
sw_im_R.py: Rosenbrock version of implicit method \
sw_trbdf2_R.py: Rosenbrock version of TR-BDF2 \
sw_im_E.py: Eisenstat-Walker version of implicit method \
sw_trbdf2_E.py: Eisenstat-Walker version of TR-BDF2 \
sw_convergence.py: calculate L2-norm between IM (small dt) and other methods to understand convergence rate as dt varies\
sw_convergence_s.py: calculate L2-norm between IM (small dt) and other methods to understand convergence rate as ref_level varies\
sw_energy.py: Analyse the energy as each time (also iterations and stepcount) \
sw_fourier.py: Perform Fourier analysis to better understand the damping effect of various methods \
sw_gamma.py: Analyse the stablility and effectiveness of TR-BDF2 with different gamma value \
sw_height: Plot the height of a point against time \
sw_iteration.py: Record the number of steps and iterations per step done by a method \
sw_spectral: Perform Hann Window on data \ 
sw_spatial.py: Look at distributions of energy as spatial discretization ref_level changes\
sw_vorticity.py: Record vorticity over the time period \


Supportive File: \
mg.py: Python file needed for multigrid method \
sw_create_mesh.py: Create mesh to be used by different method so that firedrake.assemble works properly

Redundant Files: \
sw_test.py: File containing IM and TRBDF2 methods to compare their L2-norm directly \
sw_conver.py: modified sw_convergence script to use sw_test.py 
sw_rosenbrock.py: Similar to sw_convergence but for rosenbrock method \
sw_vort.py: Similar to sw_height but for vorticity \
