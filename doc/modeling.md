Integrative modeling of ADP-actin, gelsolin and C-terminal actin-binding domain of tropomodulin {#modeling}
===============================================================================================

[TOC]

# Intrduction {#modelintro}

Here, we demonstrate integrative modeling using the PMI interface by modeling
the complex of actin and tropomodulin-gelsolin chimera using SAXS, EM,
crosslinking, crystal structures of the individual domains, and physical
principles. This complex was solved via X-ray crystallography at 2.3 Å
resolution (PDB: 4PKI). We use this structure to simulate biophysical data
and assess the accuracy of the modeled complexes. In this simple exercise,
we assume that we have a crystal structure of only the actin-gelsolin
interface and would like to find the tropomyosin-actin binding interface.
The entire modeling protocol is summarized in the four-stage diagram below:

\image html four_stage.png "The four stages of integrative modeling of actin-tropomodulin-gelsolin" width=800px

# Tutorial data files {#datafiles}

To follow this tutorial, first download the data files, either by
[cloning the GitHub repository](https://github.com/salilab/actin_tutorial)
or by [downloading the zip file](https://github.com/salilab/actin_tutorial/archive/master.zip).

# Gathering and preparing information {#gatherdata}

All data is contained in subfolders of the `data` directory of the tutorial.

## Structural data from the PDB {#gatherpdb}

The crystal structure 4PKI is used to set the atomic coordinates for each
of the domains in the FASTA sequence that determines the composition of each
biomolecule as well as the coordinates for tropomyosin and the actin-gelsolin
complex.

## Chemical Crosslinks {#gatherxlink}
Thirty-three simulated crosslinks were generated from a random subset of lysine residues whose CA-CA distances are under 25 Å.

## Electron Microscopy {#gatherem}
A simulated EM density of the entire complex was created at 20Å resolution
using %IMP. The simulated map is approximated as a Gaussian Mixture Model
(GMM).

\note Simulated EM maps can be created map can be created in %IMP using the
      following command: `simulate_density_from_pdb <file.pdb> <output.mrc> <resolution> <a/pixel>`

## SAXS {#gathersaxs}
A simulated SAXS profile of the entire 4pki.pdb complex was created using FoXS.

## Other information {#gatherother}
We also define restraints such as excluded volume and sequence connectivity to
add chemical and physical knowledge to the modeling protocol.

\image html complex.png "Actin-gelsolin-tropomyosin complex. Top: Reference crystal structure 4PKI showing actin in green, gelsolin in red and tropomyosin in blue. Bottom: Multi-scale representation and position of the system after domain shuffling and bead relaxation. Structured domains are represented by spherical beads of 1 and 10 residues. Unstructured residues from the linker between the gelsolin and tropomyosin domains are represented as gray beads." width=400px
