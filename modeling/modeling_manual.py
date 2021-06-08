'''
Script for integrative modeling of the actin/tropomodulin binding interface
'''
# Imports
from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.saxs
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.em
import IMP.pmi.dof
import ihm.cross_linkers
import IMP.atom
import IMP.saxs

# Identify data files
pdb_dir = "../data/pdb/"
fasta_dir = "../data/fasta/"
saxs_data = "./derived_data/saxs/4pki.pdb.0.15.dat"
xl_data = "./derived_data/xl/derived_xls.dat"
gmm_data = "./derived_data/em/4pki_20a_50.gmm"

# Store the FASTA sequences in a dictionary
sequences = IMP.pmi.topology.Sequences(fasta_dir + "4pkh.fasta.txt")

# Restraint weights
xl_weight = 10.0
em_weight = 10.0
saxs_weight = 10.0

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# All IMP systems start out with a Model
mdl = IMP.Model()

# In PMI, we create our top-level System within the model
sys = IMP.pmi.topology.System(mdl)

# Each System must consist of one or more States.
st = sys.create_state()

# Create the molecules within each state, where each molecule defined by
# a single peptide sequence.
# Here, we have two: Actin and a chimera of Gelsolin and Tropomyosin (geltrop)
actin = st.create_molecule("A", sequence=sequences["actin"])
geltrop = st.create_molecule("G", sequence=sequences["gelsolin-tropomyosin"])

# Add structure to the molecules from PDB files.
# You can break up structures into arbitrary domains here, as shown
# below for the two proteins in the chimera.
#
# The resulting objects (a1, a21, a22) contain a set of residues
# associated with those domains

a1 = actin.add_structure(
    pdb_dir + "4pki.pdb",
    # The chain ID in the PDB that corresponds to actin
    chain_id='A')
a21 = geltrop.add_structure(
    pdb_dir + "4pki.pdb",
    # PDB chain G is the Gelsolin-Tropomyosin chimera
    chain_id='G',
    # the residue range in the PDB that we wish to use
    res_range=(52, 177),
    # the offset between the model and PDB residue numbering
    #   For example here, residue 52 in the PDB is residue 1 in our
    #   model (based on the input FASTA sequence)
    offset=-51)
a22 = geltrop.add_structure(
    pdb_dir + "4pki.pdb",
    chain_id='G',
    res_range=(1170, 1349),
    # A discontinuity in the numbering of the PDB makes this offset a bit ugly.
    # In chain G of the PDB file, the tropomyosin domain starts at
    # residue 1170. In the actual sequence the tropomyosin domain starts
    # at residue 196.
    offset=-1170-51+196)

#####################################################
#                   REPRESENTATION                  #
#####################################################

# --------------------------
# %%%%% ADD REPRESENTATIONS TO THE SYSTEM
#
# Here, we define how we want each atom/residue in our system to be modeled.
#
# For each molecule, we choose sets of residues and add attributes to them.
#   The sets of residues can be defined by:
#       - The structured residues from PDB files defined above (a1, a21, a22)
#       - A residue range for each molecule
#         (actin[10:100] = residues 10-100 of actin)
#       - Selecting all residues with structure
#         ( geltrop.get_atomic_residues() ) or without structure
#   These sets can also be added and subtracted
#   ( geltrop[:] - geltrop.get_atomic_residues()
#     == geltrop.get_non_atomic_residues() )
#
#   Within these statements, we define the resolution(s) at which we want
#   to represent these areas and can add density components if we are
#   fitting to an EM map. We approximate the density of the model by
#   calculating Gaussian mixture models (GMMs) for each component.
#
actin.add_representation(
    a1,
    # Multi-scale representation - use beads of 1 and 10 residues per
    resolutions=[1, 10],
    # how many residues per Gaussian density particle
    density_residues_per_component=10,
    # will write a .txt and .mrc file for computed GMM (to visualize)
    density_prefix="./gmm_files/actin_gmm",
    # set to True if you want to compute the GMM each time
    density_force_compute=False,
    # Voxel size for writing the MRC file
    density_voxel_size=3.0)

actin.add_representation(
    actin[:]-a1,
    # areas without structure can only be represented at one resolution
    resolutions=[1],
    # Set up spherical gaussian densities for these particles
    setup_particles_as_densities=True)

geltrop.add_representation(
    geltrop.get_atomic_residues() and geltrop[:126], resolutions=[1, 10],
    density_residues_per_component=10,
    density_prefix="./gmm_files/gelsolin_gmm",
    density_force_compute=False,
    density_voxel_size=3.0)

geltrop.add_representation(
    geltrop.get_atomic_residues() and geltrop[144:], resolutions=[1, 10],
    density_residues_per_component=10,
    density_prefix="./gmm_files/tropomyosin_gmm",
    density_force_compute=False,
    density_voxel_size=3.0)

geltrop.add_representation(
    geltrop.get_non_atomic_residues(),
    resolutions=[1],
    setup_particles_as_densities=True)

# Build the system. Representation cannot be changed after this point!
# Any residues without a representation will be destroyed!
root_hier = sys.build()
# The Root Hierarchy of the System (root_hier) is the key node that
# we will use below to define many of our restraints.

# --------------------------
# %%%%% DEFINE THE STRUCTURAL DEGREES OF FREEDOM
#
# We need to decide what particles will move and what will stay rigid with
# respect to each other (called rigid bodies).  We also must decide the
# maximum moves step allowed for each moving particle type. These choices are
# stored in the DegreesOfFreedom object (dof)
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)

# TROPOMYOSIN RIGID BODY:
# Tropomyosin is residues 144-END in the chimera
tropo_rb_resis = geltrop[144:]
rb1 = dof.create_rigid_body(
    tropo_rb_resis,
    # the maximum translation (in angstroms) for this rigid body
    max_trans=1.0,
    # the maximum rotation of the rigid body in degrees.
    max_rot=0.5,
    # nonrigid_parts are the residues that are allowed to move within the
    # rigid body. Generally we use this for disordered loops and linkers
    # between domains. Here, we define any of the RB residues that were
    # not covered in the PDB
    nonrigid_parts=tropo_rb_resis & geltrop.get_non_atomic_residues())

# ACTIN-GELSOLIN RIGID BODY:
actgel_rb_resis = geltrop[:126]  # Gelsolin is residues 1-126 in the chimera
# We add all the actin residues (actin[:]) with the set addition operator, |=
actgel_rb_resis |= actin[:]

# Define the nonrigid residues as
nars = geltrop[:126] & geltrop.get_non_atomic_residues()
nars |= actin.get_non_atomic_residues()

rb2 = dof.create_rigid_body(
    actgel_rb_resis, max_trans=1.0, max_rot=0.5, nonrigid_parts=nars)

# Finally, we need to define movement operations for the remainder of
# the system not covered by the rigid bodies and define these as flexible
# beads with only a translation step.
fb_movers = dof.create_flexible_beads(
    geltrop.get_non_atomic_residues(), max_trans=1.0)

#####################################################
#                     RESTRAINTS                    #
#####################################################

# Restraints define functions that score the model based on
# input information.
#
# Restraint objects are first created in the definition.
# To be evaluated, the restraint object must be add_to_model().
#
# In some cases, sampled parameters for restraints must be added to the DOF
# object

# The output_objects list is used to collect all restraints
# where we want to log the output in the STAT file.
# Each restraint should be appended to this list.
output_objects = []

# -----------------------------
# %%%%% CONNECTIVITY RESTRAINT
#
# Restrains residues/particles that are collected in sequence
# This should be used for any system without an atomic force field
# (e.g. CHARMM)
# Here, we pass root_hier to apply this restraint to the entire system
cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(actin)
cr.add_to_model()           # add restraint to the model
output_objects.append(cr)   # add restraint to the output

cr2 = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(geltrop)
cr2.add_to_model()           # add restraint to the model
output_objects.append(cr2)   # add restraint to the output
# -----------------------------
# %%%%% EXCLUDED VOLUME RESTRAINT
#
# Keeps particles from occupying the same area in space.
# Here, we pass a list of both molecule chains to included_objects to
# apply this to every residue.
# We could also have passed root_hier to obtain the same behavior.
#
# resolution=1000 applies this expensive restraint to the lowest
# resolution for each particle.
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                            included_objects=[actin, geltrop],
                                            resolution=1000)
output_objects.append(evr)

# -------------------------
# %%%%% SAXS RESTRAINT
#
# Scores the SAXS data against the predicted SAXS for the input_objects
# of the model
#
# The form factor type (ff_type) for calculating the model SAXS profile
# can be either:
#   IMP.saxs.RESIDUES - Residue-level calculation - requires a model
#                       at resolution=1
#   IMP.saxs.CA_ATOMS - Residue-level calculation - requires a model
#                       at atomic resolution
#   IMP.saxs.HEAVY_ATOMS - No explicit hydrogens - requires a model
#                          at atomic resolution
#   IMP.saxs.ALL_ATOMS - Uses explicit hydrogens - requires a model
#                        at atomic resolution with hydrogens
sr = IMP.pmi.restraints.saxs.SAXSRestraint(
    input_objects=[actin, geltrop],
    saxs_datafile=saxs_data,
    weight=saxs_weight,         # scaling factor for the SAXS score
    ff_type=IMP.saxs.RESIDUES,
    # Maximum q at which to compare SAXS curves. Defaults are:
    #   ~0.15 for residue-level calculations
    #   0.4 for HEAVY_ATOMS
    #   0.5 for ALL_ATOMS (up to 1.0)
    maxq=0.15)
output_objects.append(sr)

# -------------------------
# %%%%% CROSSLINKING RESTRAINT
#
# Restrains two particles via a distance restraint based on
# an observed crosslink.
#
# First, create the crosslinking database from the input file
# The "standard keys" correspond to a crosslink csv file of the form:
#
# Protein1,Residue1,Protein2,Residue2
# A,18,G,24
# A,18,G,146
# A,50,G,146
# A,50,G,171
# A,50,G,189
#
# This restraint allows for ambiguity in the crosslinked residues,
# a confidence metric for each crosslink and multiple states.
# See the PMI documentation or the MMB book chapter for a
# full discussion of implementing crosslinking restraints.

# This first step is used to translate the crosslinking data file.
# The KeywordsConverter maps a column label from the xl data file
# to the value that PMI understands.
xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
# Here, we just use the standard keys.
xldbkc.set_standard_keys()
# One can define custom keywords using the syntax below.
# For example if the Protein1 column header is "prot_1"
# xldbkc["Protein1"]="prot_1"

# The CrossLinkDataBase translates and stores the crosslink information
# from the file "xl_data" using the KeywordsConverter.
xldb = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb.create_set_from_file(file_name=xl_data,
                          converter=xldbkc)

xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=root_hier,    # Must pass the root hierarchy to the system
    database=xldb,          # The crosslink database.
    length=25,              # The crosslinker plus side chain length
    resolution=1,           # The resolution at which to evaluate the crosslink
    slope=0.000001,         # This adds a linear term to the scoring function
                            #   to bias crosslinks towards each other
    weight=xl_weight,       # Scaling factor for the restraint score.
    linker=ihm.cross_linkers.dss)  # The linker chemistry

xlr.add_to_model()
output_objects.append(xlr)

# -------------------------
# %%%%% EM RESTRAINT
#
# Scores a model based on its cross-correlation to an EM density.
# Since cross-sorrelation is very expensive, we approximate both
# the EM map and model as a set of 3D Gaussians (done in Representation).
#
# First, collect all density particles from the model.
densities = IMP.atom.Selection(
    root_hier, representation_type=IMP.atom.DENSITIES).get_selected_particles()

emr = IMP.pmi.restraints.em.GaussianEMRestraint(
    # Evaluate the restraint using these model densities
    densities,
    # The EM map, approximated as a gaussian mixture model (GMM)
    target_fn=gmm_data,
    # a small linear restraint to pull objects towards the EM map center
    slope=0.00000001,
    # Normalizes the total density of the model wrs: EM map. Only set to true
    # if the EM map and "densities" contain the same objects.
    scale_target_to_mass=True,
    # the scaling factor for the EM score
    weight=em_weight)
emr.add_to_model()
output_objects.append(emr)

#####################################################
#                      SAMPLING                     #
#####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.

# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation=50)

# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(500)

evr.add_to_model()
emr.add_to_model()
# slx.add_to_model()
sr.add_to_model()

# Run replica exchange Monte Carlo sampling
rex = IMP.pmi.macros.ReplicaExchange0(
    mdl,
    # pass the root hierarchy
    root_hier=root_hier,
    # This allows viewing the crosslinks in Chimera
    crosslink_restraints=[xlr],
    # pass all objects to be moved ( almost always dof.get_movers() )
    monte_carlo_sample_objects=dof.get_movers(),
    # The output directory for this sampling run.
    global_output_directory='run_manual1/output/',
    # Items in output_objects write information to the stat file.
    output_objects=output_objects,
    # Number of MC steps between writing frames
    monte_carlo_steps=10,
    # set >0 to store best PDB files (but this is slow)
    number_of_best_scoring_models=0,
    # Total number of frames to run / write to the RMF file.
    number_of_frames=10000)

# Ok, now we finally do the sampling!
rex.execute_macro()

# Outputs are then analyzed using a separate analysis script.
