import numpy as np
import pandas as pd
import math
import glob
import sys
import os

sys.path.append('../lib/PMI_analysis/pyext/src/')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################

nproc = 1
top_dir = '../modeling/'
analys_dir = './'

# Check if analysis dir exists
if not os.path.isdir(analys_dir):
    os.makedirs(analys_dir)

# How are the trajectories dir names
dir_head = 'run'
out_dirs = glob.glob(top_dir+dir_head+'*/output/')

# Load module
AT = AnalysisTrajectories(out_dirs,
                          dir_name=dir_head,
                          analysis_dir = analys_dir,
                          nproc=nproc)

# Define restraints to analyze
#AT.cluster_distances =True
AT.set_analyze_XLs_restraint()
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()
AT.set_analyze_SAXS_restraint()

# Read stat files
AT.ambiguous_XLs_restraint=True
AT.read_stat_files()
AT.write_models_info()
#AT.get_psi_stats()

AT.hdbscan_clustering(['EV_sum', 'XLs_sum', 'SAXS_sum'],
                     min_cluster_size=5000)
AT.ambiguous_XLs_restraint=True
AT.summarize_XLs_info()
exit()
