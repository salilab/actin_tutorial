import numpy as np
import pandas as pd
import math
import glob
import sys
import os
import IMP

sys.path.append('../lib/PMI_analysis/pyext/src/')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################

random =0.02
top_dir = '../modeling/'


# How are the trajectories dir names
dir_head = 'run'
analys_dir = './'
out_dirs = glob.glob(top_dir+dir_head+'*/output/')

print(out_dirs)
###############################
# Extract frames
################################
# Load module
AT = AnalysisTrajectories(out_dirs,
                          dir_name=top_dir,
                          analysis_dir = analys_dir)

# Create dir
#gsms_A_dir = analys_dir+'/sample_A'
#gsms_B_dir = analys_dir+'/sample_B'

modfileA = analys_dir+'/selected_models_A_cluster-1_detailed.csv'
modfileB = analys_dir+'/selected_models_B_cluster-1_detailed.csv'

#AT.create_gsms_dir(gsms_A_dir)
#AT.create_gsms_dir(gsms_B_dir)

HA = AT.get_models_to_extract(modfileA)
HB = AT.get_models_to_extract(modfileB)

#AT.extract_models_single_rmf(HA, "selected_models_A.rmf3", top_dir, random)
#AT.extract_models_single_rmf(HB, "selected_models_B.rmf3", top_dir, random)

os.system("python ../lib/imp-sampcon/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py --sysname actin --path ../ --mode cpu_omp --cores 2 --align --density density.txt --gridsize 2.0 --gnuplot --scoreA scores_B.txt --rmfA selected_models_A.rmf3 --scoreB scores_B.txt --rmfB selected_models_B.rmf3 > clustering.out")

exit()
