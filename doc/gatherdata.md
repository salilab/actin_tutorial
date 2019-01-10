Gathering and preparing information {#gatherdata}
===================================

All data is contained in subfolders of the `data` directory of the tutorial.

# Structural data from the PDB {#gatherpdb}

The crystal structure 4PKI is used to set the atomic coordinates for each of the domains in the FASTA sequence that determines the composition of each biomolecule as well as the coordinates for tropomyosin and the actin-gelsolin complex. (Figure X, left)

# Chemical Crosslinks {#gatherxlink}
Thirty-three simulated crosslinks were generated from a random subset of lysine residues whose CA-CA distances are under 25 Å.

# Electron Microscopy {#gatherem}
A simulated EM density of the entire complex in 4pki.pdb was created at 20Å resolution using IMP. The simulated map is approximated as a Gaussian Mixture Model (GMM).

# SAXS {#gathersaxs}
A simulated SAXS profile of the entire 4pki.pdb complex was created using FoXS(Schneidman-Duhovny et al., 2010).

# Other information {#gatherother}
We also define restraints such as excluded volume and sequence connectivity to add chemical and physical knowledge to the modeling protocol. These restraints are based on XX.
