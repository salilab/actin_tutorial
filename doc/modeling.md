Integrative modeling of ADP-actin, gelsolin and C-terminal actin-binding domain of tropomodulin {#modeling}
===============================================================================================

[TOC]

# Introduction {#modelintro}

Here, we demonstrate integrative modeling using the PMI interface by modeling
the complex of actin and tropomodulin-gelsolin chimera using SAXS, EM,
crosslinking, crystal structures of the individual domains, and physical
principles. This complex was solved via X-ray crystallography at 2.3 Å
resolution (PDB: [4PKI](https://www.rcsb.org/structure/4pki)). We use this
structure to simulate biophysical data
and assess the accuracy of the modeled complexes. In this simple exercise,
we assume that we have a crystal structure of only the actin-gelsolin
interface and would like to find the tropomyosin-actin binding interface.
The entire modeling protocol is summarized in the four-stage diagram below:

\image html four_stage.png "The four stages of integrative modeling of actin-tropomodulin-gelsolin" width=800px

# Tutorial data files {#datafiles}

To follow this tutorial, first download the data files, either by
[cloning the GitHub repository](https://github.com/salilab/actin_tutorial)
or by [downloading the zip file](https://github.com/salilab/actin_tutorial/archive/main.zip).

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
using %IMP. The simulated map is approximated as a
[Gaussian Mixture Model](https://doi.org/10.1529/biophysj.108.137125) (GMM).

\note Simulated EM maps can be created map can be created in %IMP using the
      following command: `simulate_density_from_pdb <file.pdb> <output.mrc> <resolution> <a/pixel>`

## SAXS {#gathersaxs}
A simulated SAXS profile of the entire 4pki.pdb complex was created using
[FoXS](https://salilab.org/foxs/).

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

\code
|molecule_name | color | fasta_fn | fasta_id | pdb_fn | chain | residue_range | pdb_offset | bead_size | em_residues_per_gaussian | rigid_body | super_rigid_body | chain_of_super_rigid_bodies |
|actin   |green |4pki.fasta.txt|actin               |4pki.pdb|A|1,END    |0    |1|10|1|1||
|geltrop |red   |4pki.fasta.txt|gelsolin-tropomyosin|4pki.pdb|G|52,177   |-51  |1|10|1|1||
|geltrop |gray  |4pki.fasta.txt|gelsolin-tropomyosin|BEADS   |G|178,195  |-51  |1|10|1|1||
|geltrop |blue  |4pki.fasta.txt|gelsolin-tropomyosin|4pki.pdb|G|1170,1349|-1025|1|10|2|1||
\endcode

This topology file also places all domains in a single `super_rigid_body`.
This definition allows the entire complex to move as a single unit, which
is useful for fitting to the EM map.

# Constructing the modeling script {#modelscript}

The modeling script contains the entire workflow from defining the system
representation through execution of sampling. The system representation
and sampling degrees of freedom can be built manually or, as here, read
from a topology file. Restraints are added, and the sampling protocol
defined and executed.

\note The file `modeling/modeling_manual.py` contains this exact system
      built manually using PMI commands instead of a topology file.
      PMI commands allow significantly more flexibility in model design.

## Importing and building system representation {#buildsysrep}

First, we create an %IMP [Model](@ref IMP::Model) object, which stores all
components of the model. Second, we create a
[BuildSystem](@ref IMP::pmi::macros::BuildSystem) object and define the
resolutions at which residues in the structured sections will be modeled.
Here, we set resolutions of 1 and 10 residues per bead so that crosslinking
restraints can be evaluated at residue resolution and the expensive excluded
volume restraint (below) can be evaluated at the lower resolution.
Third, the topology file is read using a
[TopologyReader](@ref IMP::pmi::topology::TopologyReader) object, followed by
generating a useful list of component molecules. To this
[BuildSystem](@ref IMP::pmi::macros::BuildSystem) object, we add a state
corresponding to the representation defined in the topology file using
[bs.add_state()](@ref IMP::pmi::macros::BuildSystem::add_state).

\note To add a second state with the same topology, this line can be repeated,
      or to use a different topology, `bs.add_state(t2)` can be invoked
      with a different topology file.

\code{.py}
mdl = IMP.Model()
bs = IMP.pmi.macros.BuildSystem(mdl, resolutions=[1,10])
t = IMP.pmi.topology.TopologyReader(topology.txt)
molecules = t.get_components()
bs.add_state(t)
\endcode

We then [execute the macro](@ref IMP::pmi::macros::BuildSystem::execute_macro),
which returns the `root_hier` root hierarchy and
`dof` degrees of freedom objects, which will be used later. Within the macro,
we set the movement parameters of individual beads and rigid bodies.
Translations (`trans`) are defined in angstroms and rotations (`rot`)
in radians:

\code{.py}
root_hier, dof = bs.execute_macro(max_rb_trans=1.0,
                                  max_rb_rot=0.5, max_bead_trans=2.0,
                                  max_srb_trans=1.0, max_srb_rot=0.5)
\endcode

## Adding restraints to the model {#addrsr}

PMI contains simple interfaces for a number of %IMP restraints that model
various types of chemical and physical data and knowledge. All of these
restraints produce output, which we will collect in an `output_objects`
list. Each restraint also needs to be explicitly added to the scoring function
for sampling,using the `add_to_model()` command. We will add the restraints
to the scoring function in a specific order, discussed below. First, we
define the restraints that enforce physical and chemical principles.
(For coarse-grained models, a molecular mechanics forcefield is not
applicable. The CHARMM force field can be applied to enforce stereochemistry
on atomic models, however. See the examples in the
[IMP.atom module](@ref IMP::atom) to learn how to implement this restraint.)

The [ConnectivityRestraint](@ref IMP::pmi::restraints::stereochemistry::ConnectivityRestraint)
adds a bond between each pair of consecutive residues in each molecule.
The [ExcludedVolumeSphere](@ref IMP::pmi::restraints::stereochemistry::ExcludedVolumeSphere)
restraint is applied to the entire system and enforced at the lowest
resolution possible (indicated by resolution=1000), because this restraint
is costly to evaluate:

\code{.py}
output_objects=[]

for m in molecules:
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()
    output_objects.append(cr)

evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                     included_objects=[root_hier],
                                     resolution=1000)
\endcode

Second, we build a
[SAXSRestraint](@ref IMP::pmi::restraints::saxs::SAXSRestraint)
based on the comparison of SAXS data to the model. Since our model is
calculated at residue resolution, we calculate the SAXS profile using
residue form factors. For residue-based calculations, we compare curves out
to a q of 0.15. (Model SAXS profiles can be computed using residues, CA atoms,
heavy atoms or all atoms, depending on the resolution of the model. The
recommended `maxq` values are dependent on this choice. At residue resolution,
the fit is only valid up until q ~ 0.15; for heavy atoms q = 0.4; and for
all atoms, the fit is valid out to q = 1.0 (the maximum value).)

\code{.py}
sr = IMP.pmi.restraints.saxs.SAXSRestraint(input_objects=[root_hier],
                 saxs_datafile=saxs_data,
                 weight=0.01,
                 ff_type=IMP.saxs.RESIDUES,
                 maxq=0.15)
\endcode

To set up a crosslinking restraint, we first build a PMI
[CrossLinkDataBase](@ref IMP::pmi::io::crosslink::CrossLinkDataBase)
that uses a
[CrossLinkDataBaseKeywordsConverter](@ref IMP::pmi::io::crosslink::CrossLinkDataBaseKeywordsConverter)
to interpret a crosslink data file. At a minimum, the crosslink data file
needs four columns labeled with a key: one for each protein name and one for
each residue number of the crosslink. The standard keys are Protein1, Residue1,
Protein2, Residue2 (see `derived_xls.dat` and the `modeling.py` script for a
more in-depth explanation of crosslink keys).

\code{.py}
xl_data = "./derived_data/xl/derived_xls.dat

xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()
xldb = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb.create_set_from_file(file_name=xl_data, converter=xldbkc)
\endcode

Using this database, we can construct the crosslinking restraint. We input
the root hierarchy of the system and the database, and specify the length
of the crosslinker. The restraint can be evaluated at any resolution,
though is generally most informative at resolution = 1. The length
determines the inflection point of the scoring function sigmoid and
is generally set to 10Å+ the crosslinker length for Lys-Lys crosslinkers.

\code{.py}
xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the system root hierarchy
                CrossLinkDataBase=xldb, # The crosslink database.
                length=25,              # The crosslinker plus side chain length
                resolution=1,           # The resolution to evaluate the crosslink
                slope=0.0001,           # This adds a linear term to the score
                                        # to bias crosslinks towards each other
                weight=10)              # Scaling factor for the restraint score.

output_objects.append(xlr)
\endcode

The EM restraint is determined by calculating the overlap (cross-correlation)
between the system GMM density particles and the map GMM particles. First,
we must collect the density particles using an %IMP
[Selection](@ref IMP::atom::Selection). We then invoke the restraint using
these particles and the gmm file generated from the EM map.

\code{.py}
densities = IMP.atom.Selection(root_hier,
                 representation_type=IMP.atom.DENSITIES).get_selected_particles()

em_map = "./derived_data/em/4pki_20a_50.gmm"

emr = IMP.pmi.restraints.em.GaussianEMRestraint(
        densities,        # Evaluate the restraint using these model densities
        target_fn=em_map, # The EM map approximated as a Gaussian mixture model (GMM)
        slope=0.00000001, # a small force to pull objects towards the EM map
        scale_target_to_mass=True, # Normalizes the mass of the model wrs: EM map
        weight=100)       # the scaling factor for the EM score

output_objects.append(emr)
\endcode

## Defining the sampling protocol {#sampprot}

Sampling begins by randomizing the coordinates of the starting particles using
[shuffle_configuration](@ref IMP::pmi::tools::shuffle_configuration).
Because this randomization generally places beads of neighboring residues
far apart, we first optimize the positions of these flexible beads using
steepest descent minimization for 500 steps based on only the connectivity
restraint. We then add the balance of the scoring function terms to the
model prior to the main sampling step.

\note The shuffle algorithm fails if it cannot find a configuration without
      any overlap between components. If this happens, try increasing the
      `max_translation` parameter. Don’t set this too high as you’ll spend
      way too much time getting your system to move back together.

\code{.py}
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation=50)

dof.optimize_flexible_beads(500)

evr.add_to_model()
emr.add_to_model()
xlr.add_to_model()
sr.add_to_model()
\endcode

We implement a Monte Carlo sampling scheme with replica exchange using the PMI
[ReplicaExchange0](@ref IMP::pmi::macros::ReplicaExchange0) macro. Within
this macro, we set the directory where all output files will be placed,
`global_output_directory`, and the `number_of_frames` to generate. The final
line of the script executes the sampling macro.

\code{.py}
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
        root_hier=root_hier,           # the system root hierarchy
        crosslink_restraints= [xlr],   # This allows viewing of crosslinks in Chimera
        monte_carlo_sample_objects=dof.get_movers(), # all objects to be moved
        global_output_directory='run1/'  # Set the output directory for this run.
        output_objects=output_objects,   # Write these items to the stat file
        monte_carlo_steps=10,            # Number of MC steps between writing frames
        number_of_best_scoring_models=0, # set >0 to store best scoring PDB files
        number_of_frames=10000)          # Total number of frames to generate

rex.execute_macro()
\endcode

# Running the modeling script {#runscript}

Modeling analysis requires at least two independent sampling runs be performed.
For each run, in `modeling.py` the `global_output_directory` keyword can be
set to `run1`, `run2`, ..., `runX`.

The modeling script can be run on a single processor using the following
command:

\code{.sh}
python ../modeling.py
\endcode

or in parallel using N processors using:

\code{.sh}
mpirun -np N python ../modeling.py
\endcode

A parallel invocation of %IMP will run replica exchange with N replicas.
A serial run will run a basic Monte Carlo protocol with one replica.

Raw output will be written to the `runX/output` folder, as specified in the
replica exchange macro. Within this folder, stat files contain tabulated
statistics for each frame. In the `rmf` directory, model coordinates for
the lowest temperature replica are stored. These can be opened directly
in [UCSF Chimera](https://www.rbvi.ucsf.edu/chimera/)
and the "trajectories" observed.

Next, on to
[analysis of the actin complex models and deposition](@ref analysis).
