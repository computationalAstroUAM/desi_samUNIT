# In taurus: need to run python3 within bash

import glob
import h5py
import numpy as np
import original_read_SAGE as sage

Testing = False

outdir = '/home2/vgonzalez/outputs/desi_samUNIT/'

rootdir = '/data2/users/astevens/SAGE_output/'

##zz = 0.987
zz = 1.372

sims = ['UNITSIM1']#, 'UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
lboxes = [1000.]*len(sims) # Mpc/h
h0 = 0.6774

if Testing: sims = [sims[0]]

# Fields to read
# Type : Galaxy type (0=central, 1=satellite)
# CtreesCentralID : Central (main) halo index from consistent-trees output.
#                   This is equivalent to the MainHaloID in the ELG only SAGE h5 files.
# CentralMvir : M200c of this galaxy's main host halo [10^10 Msun/h].
# Pos : (xgal,ygal,zgal) coordinates of galaxy [cMpc/h]
# Vel : (vxgal,vygal,vzgal) velocity of galaxy [km/s]
# StellarMass : Total stellar mass of the galaxy [10^10 Msun/h]
# SfrDisk : Star formation rate of the disk (averaged over the last time-step) [Msun/yr]
# SfrBulge : Star formation rate of the bulge (averaged over the last time-step) [Msun/yr]
# BlackHoleMass : Mass of supermassive black hole [10^10 Msun/h]
# ColdGas : Cold (disc) gas mass of the galaxy [10^10 Msun / h]
# BulgeMass : Stellar mass of the bulge [10^10 Msun / h]
# DiskRadius: Exponential disk scale radius [pMpc/h]
fields=['Type','CtreesCentralID','CentralMvir','Pos','Vel',
        'StellarMass','SfrDisk','SfrBulge','BlackHoleMass',
        'ColdGas','BulgeMass','DiskRadius']
labels =[' Galaxy type (0=central, 1=satellite)',
         ' Central (main) halo index from ConsistentTrees output.',
         ' M200c of this galaxys central halo [10^10 Msun/h]',
         ' Coordinates of galaxy [cMpc/h]',
         ' Velocity of galaxy [km/s]',
         ' Total stellar mass of the galaxy [10^10 Msun/h]',
         ' Disk star formation rate (averaged over the last time-step) [Msun/yr]',
         ' Bulge star formation rate (averaged over the last time-step) [Msun/yr]',         
         ' Mass of supermassive black hole [10^10 Msun/h] ',
         ' Cold (disc) gas mass of the galaxy [10^10 Msun / h]',
         ' Stellar mass of the bulge [10^10 Msun / h]',
         ' Exponential disk scale radius in physical units [pMpc/h]']

for isim,sim in enumerate(sims):
    rootf = rootdir+sim+'/model_z'+str(zz)+'_'
    files = glob.glob(rootf+'*')
    if (len(files)<1): continue
    if Testing: files = [files[0],files[1]]
    print('{}: {} files'.format(sim,len(files)))
    
    # Prepare output file with a header
    outf = outdir+sim+'/h5_files/SAGE_z'+str(zz)+'.hdf5'
    hf = h5py.File(outf, 'w')
    head = hf.create_dataset('header',(1,))
    head.attrs[u'sim'] = sim
    head.attrs[u'box_side'] = lboxes[isim]
    head.attrs[u'units_volume'] = u'(Mpc/h)**3'
    head.attrs[u'redshift'] = zz

    # HDF5 group to store the data
    hfdat = hf.create_group('data')

    # Read all the files for a snapshot
    for ifil,fname in enumerate(files):
        for jj, field in enumerate(fields):
            alldata = sage.sageoutsingle(fname, fields=[field])
            nadd = len(alldata)
            if(nadd < 1): continue

            if (field == 'Pos' or field == 'Vel' or field == 'Spin'):
                li3 = ['x','y','z']
                for i3 in range(3):
                    data = np.zeros(shape=nadd)
                    for ii in range(nadd):
                        data[ii] = alldata[ii][0][i3]                        

                    field3 = li3[i3]+field
                    if (ifil == 0):
                        hfdat.create_dataset(field3,data=data,maxshape=(None,))
                        hfdat[field3].dims[0].label = li3[i3]+labels[jj]
                    else:
                        dset = hfdat[field3]
                        dset.resize(dset.shape[0]+nadd, axis=0)
                        dset[-nadd:,] = data
            else:
                data = np.zeros(shape=nadd)
                for ii in range(nadd):
                    data[ii] = alldata[ii][0]

                if (ifil == 0):
                    hfdat.create_dataset(field,data=data,maxshape=(None,))
                    hfdat[field].dims[0].label = labels[jj]
                else:
                    dset = hfdat[field]
                    dset.resize(dset.shape[0]+nadd, axis=0)
                    dset[-nadd:,] = data

    hf.close()
    print('Output: ',outf)
