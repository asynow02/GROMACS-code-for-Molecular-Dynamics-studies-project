Hello! 

This folder contains a .sh file with the code used for the Molecular Dynamics simulation done in open-source software GROMACS. 

References: 
Code, along with all .mdp files used, was obtained from the GROMACS tutorial 1 - Lysozyme in Water, made by Dr Justin A. Lemkul http://www.mdtutorials.com/gmx/lysozyme/index.html

Some parts were added and are not in the original code due to specific data required for analysis, such as hydrogen bonds data or Lennard-Jones and Coulombic energies data.

The only thing changed in md.mdp file was the time frame of the simulation from 1 ns to 100 ns and later for energy data, the energygrps were added.

Code in the file presents how to prepare 1ELR pdb structure - Hop TPR2A domain along with Hsp90 peptide, for the simulation and how to produce the simulation. The same code was used for other structures from the project - 2NC9, 3ESK, 2BUG etc. with small differences which are discussed within the comments. 

