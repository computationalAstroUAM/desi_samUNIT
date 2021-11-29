import glob
import numpy as np
import read_SAGE as sage

Testing = True

rootdir = '/data2/users/astevens/SAGE_output/'

zz = 0.987

sims = ['UNITSIM1']#,'UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
if Testing: sims = [sims[0]]
lboxes = [1000.]*len(sims) # Mpc/h
h0 = 0.6774

# Fields
# Type : Galaxy type (0=central, 1=satellite)
# CtreesCentralID : Central halo index from consistent-trees output.
# CentralMvir : M200c of this galaxy's central halo [10^10 Msun / h]
# Pos : (x,y,z) coordinates of galaxy [cMpc/h]
# Vel : (vx,vy,vz) velocity of galaxy [km/s]
fields=['Type','CtreesCentralID','CentralMvir','Pos','Vel']

for sim in sims:
    rootf = rootdir+sim+'/model_z'+str(zz)+'_'
    files = glob.glob(rootf+'*')
    if (len(files)<1): continue
    if Testing: files = [files[0]]

    for fname in files:       
        data = sage.sageoutsingle(fname, fields)

        #here: Write down info into a hdf5 file that appends the data from the different files, in same units as for ELGs.  
