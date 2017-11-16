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
import IMP.atom
import IMP.saxs

# Identify data files
pdb_dir = "../data/pdb/"
fasta_dir = "../data/fasta/"
saxs_data = "./derived_data/saxs/4pki.pdb.0.15.dat"
xl_data = "./derived_data/xl/derived_xls.dat"
gmm_data = "./derived_data/em/4pki_20a_50.gmm"

sequences = IMP.pmi.topology.Sequences(fasta_dir + "4pkh.fasta.txt")

# Restraint weights
xl_weight = 10.0
em_weight = 10.0
saxs_weight = 10.0


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mdl = IMP.Model()
sys = IMP.pmi.topology.System(mdl)
st = sys.create_state()

actin = st.create_molecule("A", sequence=sequences["actin"])
geltrop = st.create_molecule("G", sequence=sequences["gelsolin-tropomyosin"])

# Read in PDB files. 
a1 = actin.add_structure(pdb_dir + "4pki.pdb",
                        chain_id='A')
a21 = geltrop.add_structure(pdb_dir + "4pki.pdb",
                        chain_id='G', res_range=(52,177), offset=-51)
a22 = geltrop.add_structure(pdb_dir + "4pki.pdb",
                        chain_id='G', res_range=(1170,1349), offset=196-1170-51)


#####################################################
################### REPRESENTATION ##################
#####################################################
actin.add_representation(a1, resolutions=[1,10],
                        density_residues_per_component=10, #how much to coarsen this representation
                        density_prefix="./gmm_files/actin_gmm",         # will write a .txt and .mrc file forcomponent
                        density_force_compute=False,       # set True if you want to overwrite
                        density_voxel_size=3.0)            # set to 0 if you don't care about writing the map)

actin.add_representation(actin[:]-a1, resolutions=[1],
                        setup_particles_as_densities=True)

geltrop.add_representation(geltrop.get_atomic_residues() and geltrop[:126], resolutions=[1,10],
                        density_residues_per_component=10,
                        density_prefix="./gmm_files/gelsolin_gmm",      
                        density_force_compute=False,      
                        density_voxel_size=3.0)    

geltrop.add_representation(geltrop.get_atomic_residues() and geltrop[144:], resolutions=[1,10],
                        density_residues_per_component=10,
                        density_prefix="./gmm_files/tropomyosin_gmm",      
                        density_force_compute=False,      
                        density_voxel_size=3.0)          

geltrop.add_representation(geltrop.get_non_atomic_residues(), resolutions=[1],
                            setup_particles_as_densities=True)

# Build the system. Representation cannot be changed after this point!
root_hier = sys.build()

# Create Rigid Bodies
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
tropo_rb_resis = geltrop[144:] 
actgel_rb_resis = geltrop[:126]
actgel_rb_resis |= actin[:]

#print geltrop.get_atomic_residues()

#exit()

for i in IMP.pmi.tools.input_adaptor(actgel_rb_resis,
                                                'all',
                                                flatten=True):
    #print i, i.get_particle()
    j=1

#print actgel_rb_resis
nars = geltrop[:126] & geltrop.get_non_atomic_residues()
nars |= actin.get_non_atomic_residues()
rb1 = dof.create_rigid_body(tropo_rb_resis, max_trans=1.0, max_rot=0.5, nonrigid_parts=tropo_rb_resis & geltrop.get_non_atomic_residues())
rb2 = dof.create_rigid_body(actgel_rb_resis, max_trans=1.0, max_rot=0.5, nonrigid_parts=nars)
fb_movers = dof.create_flexible_beads(geltrop.get_non_atomic_residues(),max_trans=1.0)

#####################################################
##################### RESTRAINTS ####################
#####################################################
output_objects = []
# -------------------------
# Connectivity Restraint
cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(geltrop)
cr.add_to_model()
output_objects.append(cr)

# -------------------------
# Excluded Volume Restraint
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=[actin, geltrop])
evr.add_to_model()
output_objects.append(evr)

# -------------------------
# SAXS Restraint
sr = IMP.pmi.restraints.saxs.SAXSRestraint(input_objects=[actin, geltrop], 
                                            saxs_datafile=saxs_data,
                                            weight=saxs_weight,
                                            ff_type=IMP.saxs.RESIDUES)
sr.add_to_model()
output_objects.append(sr)

# -------------------------
# Crosslinking Restraint

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
# Other keys can be added, see the PMI documentation for a full list of keywords

xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()
xldb = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb.create_set_from_file(file_name=xl_data,
                            converter=xldbkc)

xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                    root_hier=root_hier,
                                    CrossLinkDataBase=xldb,
                                    length=25,      # The crosslinker plus side chain length
                                    resolution=1,   # The resolution at which to evaluate the crosslink
                                    slope=0.0001,     # This adds a linear term to the scoring function 
                                                    # to bias crosslinks towards each other
                                    weight=xl_weight)

xlr.add_to_model()
output_objects.append(xlr)

xlr.evaluate()

# -------------------------
# EM Restraint
densities = IMP.atom.Selection(root_hier,representation_type=IMP.atom.DENSITIES).get_selected_particles()
emr = IMP.pmi.restraints.em.GaussianEMRestraint(
     densities,
     target_fn=gmm_data,  # created by user, see top of file
     slope=0.00000001,                      # a small number, helps drag bits into map
     scale_target_to_mass=True,           # if the model is the same size as map, usually set to True
     weight=em_weight)                         # the data weight
emr.add_to_model()
output_objects.append(emr)

#####################################################
###################### SAMPLING #####################
#####################################################
for i in IMP.atom.Selection(root_hier).get_selected_particles():
    print(i, IMP.core.RigidBody.get_is_setup(i), rb1, rb2)
# First shuffle all particles
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation=30)


# Quickly move all flexible beads into place
#dof.optimize_flexible_beads(100)
print(xlr.evaluate())
print(emr.evaluate())
print(sr.evaluate())

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=root_hier,                          # pass the root hierarchy
                                    crosslink_restraints=[cr, xlr],               # This allows viewing the crosslinks in Chimera
                                    monte_carlo_sample_objects=dof.get_movers(),  # pass all objects to be moved ( almost always dof.get_movers() )
                                    global_output_directory='output/',
                                    output_objects=output_objects,          # Items in output_objects write information to the stat file.
                                    monte_carlo_steps=10,                   # Number of MC steps between writing frames
                                    number_of_best_scoring_models=0,        # set >0 to store best PDB files (but this is slow)
                                    number_of_frames=100)                   # Total number of frames to run / write to the RMF file.

# Ok, now we finally do the sampling!
rex.execute_macro()

print(xlr.evaluate())
print(emr.evaluate())
print(sr.evaluate())
