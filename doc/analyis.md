Analysis of the actin complex models and deposition {#analysis}
===================================================

[TOC]

Analysis is performed using scripts located in `analysis/scripts`.
The already-generated sampling output will be analyzed here; it is contained
in the folders `modeling/run1` and `modeling/run2`.

Analysis is performed in a new directory: `analysis/tutorial_analysis`.

# Filtering good scoring models {#filtering}

The `select_good_scoring_models.py` script filters models based on score
and parameter thresholds. In this script, required flags are: `-rd`, which
specifies the directory containing sampling output folders; `-rp`, which
defines the prefix for the sampling output folders; `-sl`, which defines
the stat file keywords (see below) that we wish to filter on; `-pl`,
which specifies the keywords that will be written to the output file;
`-alt` and `-aut`, which specify, respectively, the lower and upper threshold
for each keyword in `-slt` hat define acceptance. The `-mlt`
and `-mut` keywords, which are optional, define thresholds for restraints made
of multiple components (such as crosslinks).

\note A list of acceptable stat file keywords can be determined by running
      `../scripts/plot_stat.py ./path/to/stat/file -pk`

Here, we first use crosslink satisfaction as an initial filtering criterion
because we usually have an *a priori* estimate of the false positive rate
and/or cutoff distance (for scores whose thresholds are not known *a priori*,
one can perform a multi-stage filtering process as outlined in the above
protocol). For this simulated system, we only accept models with 100%
satisfaction of crosslinks by setting both `-alt` and `-aut` to 1.0.
A crosslink is satisfied if the distance is between 0.0 and 30.0 Å, as
delineated by the `-mlt` and `-mut` keywords, respectively. We specify
that connectivity, crosslink data score, excluded volume, EM, SAXS and
total scores be printed as well.

\code{.sh}
python ../scripts/select_good_scoring_models.py -rd ../../modeling -rp run -sl "CrossLinkingMassSpectrometryRestraint_Distance_" -pl ConnectivityRestraint_None CrossLinkingMassSpectrometryRestraint_Data_Score ExcludedVolumeSphere_None GaussianEMRestraint_None SAXSRestraint_Score Total_Score -alt 1.0 -aut 1.0 -mlt 0.0 -mut 30.0
\endcode
