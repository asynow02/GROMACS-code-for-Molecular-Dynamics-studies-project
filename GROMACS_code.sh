
#For some of the structures, like 1ELR or 3ESK, the first step was to remove the water "HOH" molecules
#that was simply done by using the command below and saving the new structure in a seperate PDB file

#On top of that both 1ELR and 3ESK had NI ions, but they were manually removed from the pdb files

grep -v HOH 1ELR.pdb > 1ELR_clean.pdb 

#After the structure was cleaned up, the following command was used to generate topology file, a position restraint file and post processed structure .gro file
#the water model was specified, in this case spce. Additionally, -ter option was used to select manually
#N terminus of the peptide (in this case 'None' option). Since the peptide was capped, that had to be done because GROMACS automatically 
#attempts to construct a free amino group onto acetyl. That option though, was only used for structures which had Hsp90 or Hsc70 peptides

gmx pdb_2gmx -f 1ELR_clean.pdb -o 1ELR_processed.gro -water spce -ter 

#Additionally, for 2NC9 and 2BUG, which are NMR structures, option -ignh was used too, to ignore H atoms 

#After processing the command above, GROMACS asks to choose force field, in this case option 15: OPLS-AA was used

#Next step was to define simulation box, post processed structure file was used as an input file,
# -c option centered the protein in the box, -d 1.0 option put it at least 1.0 nm from the edge of the box,
# - bt option specifies the type of the box which in this case was cubic 

gmx editconf -f 1ELR_processed.gro -o 1ELR_newbox.gro -c -d 1.0 -bt cubic 

#The box created was then filled up with water molecules, by using -cs - GROMACS solvent configuration 
#and spc216.gro which is generic equilibrated 3-point solvent model box.
#Output is saved into new .gro file and the topology is updated 

gmx solvate -cp 1ELR_newbox.gro -cs spc216.gro -o 1ELR_solv.gro -p topol.top 

#Below command generates a .tpr file needed for adding ions to the structure. Some structures 
#have either positive or negative net charge. Adding ions result in making the structure neutral which is needed for obtaining the simulation 

gmx grompp -f ions.mdp -c 1ELR_solv.gro -p topol.top -o ions.tpr 

# -pname and -nname options define the positive and negative ions which are going to be added to the structure

gmx genion -s ions.tpr -o 1ELR_solv_ions.gro -p topol.top -pname NA -nname CL -neutral 

#After processing command above, GROMACS asked to choose a group for embedding ions - 13 SOL was selected to avoid 
#replacing the parts of the protein with ions

#Next step was to perform Energy Minimization, which ensures the system is prepared for simulation

gmx grompp -f minim.mdp -c 1ELR_solv_ions.gro -p topol.top -o em.tpr

#gmx mdrun was used to carry out any operations, such as Energy Minimization, or the actual simulations
#Four different files were generated - em.log, em.edr, em.trr and em.gro 

gmx mdrun -deffnm em 

#gmx energy is used to produce .xvg files, which is helpful in vizualization of different steps and ensuring everything was done correctly
#in this case option 10 0 were chosen for Potential 

gmx energy -f em.edr -o potential_1ELR.xvg 

#The system was then equilibrated under NVT ensemble to ensure the temperature is stabilized 

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr 

gmx mdrun -deffnm nvt 

gmx energy -f nvt.edr -o temperature_1ELR.xvg 

#Option 16 0 was chosen for temperature results 

#The system was then equilibrated again under NPT ensemble to ensure the temperature, pressure and number of particles are all constant

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr 

gmx mdrun -deffnm npt

gmx energy -f npt.edr -o pressure_1ELR.xvg 

#Option 18 0 was chosen for presure results

gmx energy -f npt.edr -o density_1ELR.xvg 

#Option 24 0 was chosen for density results 

#The MD simulation is then produced and saved 

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_2_1.tpr 

gmx mdrun -deffnm md_0_2_1

#The simulation file is then processed to account any priodicity in the system by using gmx trjconv 

gmx trjconv -s md_0_2_1.tpr -f md_0_2_1.xtc -o md_0_2_1_noPBC.xtc -pbc mol -center 

#Option 1 Protein is used as the group to be centered and 0 System for output 

#RMSD is then calculated by using gmx rms function, option -tu ns specifies that the time will be in nanosecond unit

gmx rms -s md_0_2_1.tpr -f md_0_2_1_noPBC.xtc -o rmsd_1ELR.xvg -tu ns 

#Same is used for RMSD in relation to the crystal structure 

gmx rms -s em.tpr -f md_0_2_1_noPBC.xtc -o rmsd_xtal_1ELR.xvg -tu ns

#In both cases option 4 Backbone is used for the least-swaures fit and the RMSD calculation group

#The radius of gyration is calculated by using gmx gyrate function 

gmx gyrate -s md_0_2_1.tpr -f md_0_2_1_noPBC.xtc -o gyrate_1ELR.xvg

#Option 1 protein was chosen 

#Following steps were done in order to obtain energy and hydrogen bonds information from the simulation file
#First the index .ndx files were created by using gmx make_ndx function.

gmx make_ndx -f em.gro -o index_carboxylate_clamp.ndx 

gmx make_ndx -f em.gro -o index_peptide.ndx 

gmx make_ndx -f em.gro -o index_rest.ndx 

#Here three different options for index file creation are shown as an example, that was obviously different 
#depending if structure included peptide or not
#As an example for the index_peptide.ndx file, 'a 1089-1131' was used, which are the numbers of the atoms from the Hsp90 peptide present in the structure
#The same was done for the carboxylate clamp and the rest of TPR domain residues, the atom numbers were provided 
#Depending on the structure that numbers will be different so it is always good to check in the em.gro file 
#later the labels were changed for easier identification by using 'name' and to exit 'q'

#indexes were then joined together using cat 

cat index_carboxylate_clamp.ndx index_peptide.ndx > index_peptide_clamp.ndx 

cat index_carboxylate_clamp.ndx index_rest.ndx > index_rest_clamp.ndx 

#hydrogen bonds data was generated using the simulation file, the appropriate index (for example the index containing information about atoms from the carboxylate
#clamp residues and the rest of the TPR domain, was used in the function to analyse the bonds between them)

gmx hbond -f md_0_2_1.xtc -s md.tpr -n index_rest_clamp.ndx -num hbond_rest_clamp_1ELR.xvg 

gmx hbond -f md_0_2_1.xtc -s md.tpr -n index_peptide_clamp.ndx -num hbond_peptide_clamp_1ELR.xvg

#to analyse the Coulombic and Lennard-Jones energies between atoms from specific residues, the simulations had to be rerun, after adding indexes and
#adding to the md.mdp energygrps line which included all the names of the groups or residues from the appropriate index

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index_rest_clamp.ndx -o md_0_2_clamp_rest.tpr 

gmx mdrun -deffnm md_0_2_clamp_rest -rerun md_0_2_1.xtc 

#here simply the energies were checked between the clamp residues and the rest of the TPR domain, but that was also done for individual carboxylate
#clamp residues and the peptide for example, to do that the appropriate index had to be generated, the md.mdp file energygrps specified 
#and then correct option in gmx energy selected - for the example below options for Coul-SR and LJ-SR between the carboxylate clamp and the rest of the TPR domain residues were selected

gmx energy -f md_0_2_clamp_rest.edr -o interaction_energy_clamp_rest_1ELR.xvg 

#additionally, RMSF was calculated using command below 

gmx rmsf -f md_0_2_1.xtc -s md.tpr -n index_TPR_domain.ndx -o rmsf_1ELR.xvg -res 

#index option was used, since rmsf of only TPR domain was needed for comparison, the index was produced 
#the same way as showed before, but along with choosing only atoms from TPR domain, 'a 1-1088 & a CA' was used to choose only Ca backbone atoms

#the new option created before 'a 1-1088 $ a CA', was chosen after running the command














