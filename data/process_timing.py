# Script to read in a list of processing times and get the average.

import numpy as np

atom_times_kd = np.loadtxt("./atom_timing_500cm_nsh.txt")
atom_times_hash = np.loadtxt("./atom_timing_500cm_nsh_hashing.txt")
print "Average insertion time for AtomMap with kdtree was %f seconds." % atom_times_kd.mean()
print "Average insertion time for AtomMap with hashing was %f seconds." % atom_times_hash.mean()
