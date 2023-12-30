# amyloid_beta_Kuznetsov

Introduction
This supplement provides the MATLAB source code accompanying the manuscript 
"Computational investigation of the effect of reduced dynein velocity and reduced cargo diffusivity on slow axonal transport" 
by Ivan A. Kuznetsov and Andrey V. Kuznetsov

Last updated: January 2023

Please put the files in the same folder, one code produces data output (.mat file) that the other code reads,

Finding values of parameters in Table S3 that give the best fit with experimental results for MAP1B protein transport. 
The global minimum was found by using routine MULTISTART with a local solver FMINCON, which can be found in Matlab’s 
Optimization Toolbox. 5000 randomly selected starting points in the parameter space were used to start the optimization 
procedure.

test_MultiStart_hpc76.m funerr5_6_hpc76.m


Numerical solution of full SAT model given by Eqs. (1)-(5) with boundary conditions (15a,b) and (16a,b).
This code is used to generate Figs. 3a, 3b, 4a, 4b in the paper and
Figs. S1a, S1b, and S1c in the Supplementary Materials.

fig5_6_main.m funerr5_6_figure.m


Numerical solution of the perturbation equations (S3)-(S7) for vr->0 with boundary conditions (S8a,b) and (S9). 
This code is used to generate Figs. 5a and 5b in the paper and Figs. S2a and S2b in the Supplementary Materials.

fig5_6_main_vr_0.m funerr5_6_figure_vr_0.m



Numerical solution of the perturbation equations (S10)-(S14) for Dfree->0 with boundary conditions (S15) and (S16). 
This code is used to generate Figs. 6a and 6b in the paper and
Figs. S3a and S3b in the Supplementary Materials.

fig5_6_main_Dfree_0.m funerr5_6_figure_Dfree_0.m



Analytical solution of the perturbation equations (S17)-(S21) for vr->0 and Dfree->0 with boundary condition (S22).
This code is used to generate Figs. 7a, 7b, and 8 in the paper.

fig5_6_main_vr_0_Dfree_0.m funerr5_6_figure_vr_0_Dfree_0.m
