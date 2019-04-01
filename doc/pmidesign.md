PMI software design {#pmidesign}
===================

[TOC]

# Goal {#goal}

The goal of the Python Modeling Interface (PMI)
is to allow structural biologists with limited programming
expertise to determine the structures of large protein complexes
by following the [IMP four-step integrative modeling protocol](@ref procedure).
PMI is a top-down modeling system that relies on a series of macros and
classes to simplify encoding of the modeling protocol, including designing the
system representation, specifying scoring function, sampling alternative
structures, analyzing the results, facilitating the creation of
publication-ready figures, and depositing [into PDB-Dev](@ref deposition).
PMI exchanges the high flexibility of %IMP for ease-of-use, all within one short
Python script (<100 lines). Despite its simplicity in creating standard
modeling workflows, PMI is powerful and extensible - it is built on %IMP and
creates native %IMP objects, which means that the advanced user can customize
many aspects of the modeling protocol. Below, we outline each stage of the
modeling process as performed in PMI.

\image html pmi_hierarchy.png "PMI Hierarchy. PMI is constructed as a top-down hierarchy beginning with a System. A system can contain one or more states with each state being a different conformation, composition or time-ordered step in the system. Each state is comprised of one or more molecules that may have one or more copies per molecule. At the final level, each molecule is represented at one or more resolutions." width=600px

# Gathering information {#gatherinfo}

Information about a system that we wish to model includes everything that we
directly observe, can infer through comparison to other systems, and
fundamental physical principles. Experimental data that are commonly utilized
in integrative modeling include X-ray crystal structures, EM density maps,
NMR data, chemical crosslinks, yeast two-hybrid data,and Förster resonance
energy transfer (FRET) measurements. Atomic resolution information may be
applied directly as structural restraints from atomic statistical potentials
and molecular mechanics force fields or derived from comparative modeling
programs such as [MODELLER](https://salilab.org/modeller/) and
[PHYRE2](https://doi.org/10.1038/nprot.2015.053). Each piece of information
can be utilized within the modeling procedure in one or more of five distinct
ways: defining model **representation**, defining **sampling** space/degrees of
freedom, **scoring** models during sampling, **filtering** models post-sampling,
and **validating** completed models.

# System and data representation {#systemdatarep}

The representation of the system defines the structural degrees of freedom
that will be sampled and is designed based on the information at hand. We can
utilize a multi-scale representation, where model components can be modeled
at one or more different resolutions commensurate with the information content
at that site. For example, a domain described by a crystal structure can be
represented at atomic resolution and a disordered segment can be represented
as a string of spherical beads of 10 residues each. In addition, non-particle
based representations, such as Gaussian mixture models (GMMs),
can also be used; for example, in EM density fitting. Choosing a
representation reflects a compromise between the need for details required
by the biological application of the model and the need for coarseness
required by limited computing power.

\image html representation.png "Ways of representing a single biomolecule. A: The complexity of a molecular system can be represented in four ways. An ensemble of states describes the structural heterogeneity around a single solution. Multiplestates are used to describe systems that exist in multiple thermodynamic wells. The system can be modeled at a multitude of scales commensurate with the different types of information known about it. Finally, individual states can be time-ordered, allowing for the modeling of the transition rates between them. B: Multiple representations can be simultaneously applied to the same biomolecule so that information of various types can be applied at the proper scale and form. The molecule is first defined by its sequence connectivity. Flexible beads comprising one or more residues are commonly applied to loops where no high-resolution structure is available. Areas that have high resolution structure can be modeled by spherical beads of 1 residue for the evaluation of residue-specific information such as chemical crosslinks or NMR distance restraints. 10-residue beads are generally used to model lower resolution information such as SAXS data and the excluded volume restraint. The molecule can be represented as a Gaussian mixture model for comparisons to EM densities. IMP and PMI can utilize all of these representations simultaneously in a multi-scale model." width=600px
