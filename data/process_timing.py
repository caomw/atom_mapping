# Script to read in a list of processing times and get the average.

import numpy as np

atom_times = np.loadtxt("./revision/lbl_atom_1m_timing_sdf.txt")
#octo_times = np.loadtxt("./revision/nsh_octo_1m_timing.txt")
print "Average insertion time for AtomMap was %f seconds." % atom_times.mean()
#print "Average insertion time for OctoMap was %f seconds." % octo_times.mean()
