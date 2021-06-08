'''
Script for integrative modeling of the actin/tropomodulin binding interface
using a topology file
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
saxs_data = "./derived_data/saxs/4pki.pdb.0.15.dat"
xl_data = "./derived_data/xl/derived_xls.dat"
gmm_data = "./derived_data/em/4pki_20a_50.gmm"

# Restraint weights
xl_weight = 10.0
em_weight = 1000.0
saxs_weight = 0.01

# Topology File
topology_file = "topology.txt"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# All IMP systems start out with a Model
mdl = IMP.Model()

# Read the topology file for a given state
t = IMP.pmi.topology.TopologyReader(topology_file)

# Create a BuildSystem macro to and add a state from a topology file
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(t)

# executing the macro will return the root hierarchy and degrees
# of freedom (dof) objects
root_hier, dof = bs.execute_macro()

# It's useful to have a list of the molecules.
molecules = t.get_components()


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
# (e.g. CHARMM). We apply the restraint to each molecule.

for m in root_hier.get_children()[0].get_children():
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()
    output_objects.append(cr)

# -----------------------------
# %%%%% EXCLUDED VOLUME RESTRAINT
#
# Keeps particles from occupying the same area in space.
# Here, we pass a list of both molecule chains to included_objects to
# apply this to every residue.
# We could also have passed root_hier to obtain the same behavior.
#
# resolution=1000 applies this expensive restraint to the lowest resolution
# for each particle.
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
    included_objects=[root_hier], resolution=1000)
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
    input_objects=[root_hier],
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
    slope=0.0001,           # This adds a linear term to the scoring function
                            #   to bias crosslinks towards each other
    weight=xl_weight,       # Scaling factor for the restraint score.
    linker=ihm.cross_linkers.dss)  # The linker chemistry


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
    #   if the EM map and "densities" contain the same objects.
    scale_target_to_mass=True,
    # the scaling factor for the EM score
    weight=em_weight)

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

# Now, add all of the other restraints to the scoring function to
# start sampling
evr.add_to_model()
emr.add_to_model()
xlr.add_to_model()
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
    global_output_directory='run2/output/',
    # Items in output_objects write information to the stat file.
    output_objects=output_objects,
    # Number of MC steps between writing frames
    monte_carlo_steps=10,
    # set >0 to store best PDB files (but this is slow)
    number_of_best_scoring_models=0,
    # Total number of frames to run / write to the RMF file.
    number_of_frames=50000)

# Ok, now we finally do the sampling!
rex.execute_macro()

print(sr.evaluate())
print(emr.evaluate())
print(xlr.evaluate())


# Outputs are then analyzed in a separate analysis script.
