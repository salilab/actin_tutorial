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

# Defining system representation and degrees of freedom in the topology file {#topology}

The model representation (e.g., bead size and rigid bodies) can be set within
the topology file. The topology file is a pipe-delimited format with each
line specifying a separate domain and keyword values determining how the
domain is represented. For more information on the topology file and its
contents, see [the PMI documentation](@ref IMP::pmi::topology::TopologyReader).

The topology file for this tutorial, shown below, is found at
`modeling/topology.txt`. Here, the system is subdivided into four distinct
domains: one each for the three structured domains (actin, gelsolin, and
tropomyosin) and one consisting of the 18-residue engineered linker between
gelsolin and tropomyosin. The first domain, the entire actin molecule, is
colored green and contains the entirety of `chain` A from 4pki.pdb.
A `bead_size` of 1 residue per bead is assigned to any unmodeled section,
i.e. not present in the PDB file (spherical beads are applied to every 10
residues with smaller beads applied to loops of smaller length). A GMM is
approximated using 10 residues per Gaussian. This domain is assigned to
`rigid_body` 1. The second domain, the gelsolin portion of the chimera, is
constructed by selecting the `residue_range` 52-177 of `chain` G.
These residues, however, are numbered 1-126 in the FASTA file, therefore a
`pdb_offset` of -51 must be added. This domain is also assigned to `rigid_body`
1 to preserve the actin/gelsolin interface. The third domain is the linker,
whose residues have no structure associated with them; thus, they are given a
`pdb_fn` of BEADS with a `bead_size` of 1 (these residues are also assigned to
`rigid_body` 1 to improve sampling; all beads within rigid bodies are, by
default, allowed to be flexible). The final domain, tropomyosin, is built
similarly to gelsolin and assigned to `rigid_body` 2,
since we would like to sample its position separate of the rest of the complex.

@code
|molecule_name | color | fasta_fn | fasta_id | pdb_fn | chain | residue_range | pdb_offset | bead_size | em_residues_per_gaussian | rigid_body | super_rigid_body | chain_of_super_rigid_bodies |
|actin   |green |4pki.fasta.txt|actin               |4pki.pdb|A|1,END    |0    |1|10|1|1||
|geltrop |red   |4pki.fasta.txt|gelsolin-tropomyosin|4pki.pdb|G|52,177   |-51  |1|10|1|1||
|geltrop |gray  |4pki.fasta.txt|gelsolin-tropomyosin|BEADS   |G|178,195  |-51  |1|10|1|1||
|geltrop |blue  |4pki.fasta.txt|gelsolin-tropomyosin|4pki.pdb|G|1170,1349|-1025|1|10|2|1||
@endcode
