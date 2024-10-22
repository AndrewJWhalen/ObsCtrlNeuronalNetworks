# ObsCtrlNeuronalNetworks
Code archive of Matlab and Mathematica scripts required to replicate numerical results published in PRX: "Observability and Controllability of Nonlinear Networks: The Role of Symmetry" (DOI: 10.1103/PhysRevX.5.011005)

This is a code archive for the paper: 

Observability and Controllability of Nonlinear Networks: The Role of Symmetry.
Andrew J. Whalen, Sean N. Brennan, Timothy D. Sauer, and Steven J. Schiff
Physical Review X, 2014

version: 11/18/2014
Andrew Whalen
email: andrew.whalen@yale.edu (Updated email as of 1/1/2023)

Unzip the file in a new directory folder.

1) Running "PlotResults.m" in Matlab will re-create Figure 2 from the paper from the results output data included in the "MathematicaData" folder.

2) Inside the "MatlabData" folder "FN3NodeSimulationData.m" is a Matlab script that will enable the creation of the various simulated data trajetories of the 3-node motifs in the paper.

3) "CtrlbMotif1.m" and "ObsMotif1.m" are Mathematica batch run files which will load the simulated data contained in the "MatlabData" folder (created from step 2) above) and compute the symbolic Controllability and Observability matrices respectively and numerically evaluate these matrices along each point of the simulated trajectories. These scripts are optimized for running on a supercomputing cluster.

4) Any output results from step 3) above may be plotted using the script in step 1) above with the proper selection of the control parameters at the beginning of the code. 

We offer these codes so that others may replicate and build upon our findings. Please cite the PRX paper (DOI: 10.1103/PhysRevX.5.011005) from which these original results were published if you make use of these codes, or modify them for your own data and research.

Andrew Whalen,
Center for Neural Engineering,
The Pennsylvania State University,
November, 2014
