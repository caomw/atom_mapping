# Script to read in a list of processing times and get the average.

import numpy as np

octo_times = np.loadtxt("./octomap_timing_810cm_lbl.txt")
atom_times = np.loadtxt("./atom_timing_500cm_lbl.txt")
print "Average insertion time for OctoMap was %f seconds." % octo_times.mean()
print "Average insertion time for AtomMap was %f seconds." % atom_times.mean()
