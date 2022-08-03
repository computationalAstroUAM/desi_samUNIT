# Create a hdf5 file with the selection of galaxies and properties 
import os.path, sys   
import h5py
import numpy as np
from Cosmology import set_cosmology,H
from iotools import check_file

Testing = False

h5file = False

zz = 0.987

sims = ['UNITSIM1']#,'UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
lboxes = [1000.]*len(sims) # Mpc/h                                                          
h0 = 0.6774

############################################
unitdir = '/data8/aknebe/UNITSIM/ELGs_DESI/'
outdir = '/home2/vgonzalez/out/desi_samUNIT/'
############################################

props = ['mass','sfr','lo2']
labelp = ['log10(M/Msun/h)','log10(SFR/Msun/Gyr)', 'log10(L[OII]/h^-2 erg/s)']

if Testing: sims = [sims[0]] ; props = [props[0]]

redshift = str(zz).replace('.','_')

# Get Hubble constant at zz
set_cosmology(omega0=0.3089,omegab=0.02234/h0/h0,h0=h0,
              universe="Flat",include_radiation=False)
Hz = H(zz)

for isim,sim in enumerate(sims):
    lbox = lboxes[isim]

    # Prep output directories
    if h5file:
        outpath = outdir+sim+'/h5_files/'
    else:
        outpath = outdir+sim+'/ascii_files/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        
    # Read the data from the sim
    ff = unitdir+sim+'/'+sim+'_model_z'+str(zz)+'_ELGs.h5'
    if (not check_file(ff)):  continue
    f = h5py.File(ff,'r')

    ifirst=0 ; ilast=f['Mstar'].shape[0]
    if Testing:
        ilast = 50000

    mhalo = f['Mhalo'][ifirst:ilast] # Msun/h
    xgal = f['Xpos'][ifirst:ilast] # Mpc/h
    ygal = f['Ypos'][ifirst:ilast] # Mpc/h
    zgal = f['Zpos'][ifirst:ilast] # Mpc/h
    xgalz = f['Xpos'][ifirst:ilast] + f['Xvel'][ifirst:ilast]*(1+zz)/Hz
    ygalz = f['Ypos'][ifirst:ilast] + f['Yvel'][ifirst:ilast]*(1+zz)/Hz
    zgalz = f['Zpos'][ifirst:ilast] + f['Zvel'][ifirst:ilast]*(1+zz)/Hz
    mass = f['Mstar'][ifirst:ilast] # Msun/h
    sfr = f['SFR'][ifirst:ilast]*10**9 # Msun/h/Gyr
    lo2 = 10**f['logLOII_3727_att'][ifirst:ilast] +\
          10**f['logLOII_3729_att'][ifirst:ilast]

    # Correct for periodic boundary effects
    ind = np.where(xgalz>lbox) ; xgalz[ind] = xgalz[ind] - lbox
    ind = np.where(xgalz<0.)   ; xgalz[ind] = xgalz[ind] + lbox
    ind = np.where(ygalz>lbox) ; ygalz[ind] = ygalz[ind] - lbox
    ind = np.where(ygalz<0.)   ; ygalz[ind] = ygalz[ind] + lbox
    ind = np.where(zgalz>lbox) ; zgalz[ind] = zgalz[ind] - lbox
    ind = np.where(zgalz<0.)   ; zgalz[ind] = zgalz[ind] + lbox

    # Make selections for different cuts
    for ip,prop in enumerate(props):
        vals = np.full(len(mhalo),-999.)
        if (prop == 'mass'): #here vals in adequate units
            array = mass
            ind = np.where((mhalo>0) & (array>0))
            vals[ind] = np.log10(array[ind])
            array = []
        elif (prop =='sfr'):
            array = sfr
            ind = np.where((mhalo>0) & (array>0))
            vals[ind] = np.log10(array[ind])
            array = []
        elif (prop == 'lo2'):
            array = lo2
            ind = np.where((mhalo>0) & (array>0))
            vals[ind] = np.log10(array[ind]) + 2*np.log10(h0) # erg/s/h**2 
            array = []

        # Read the prop cuts
        ndfile = outdir+sim+'/ngal_'+prop+'_cuts_z'+redshift+'.dat'
        nds, cuts = np.loadtxt(ndfile, unpack=True)
        if np.isscalar(nds):                                                                          
            nds = np.array([nds]) ; cuts = np.array([cuts])   

        for ic,cut in enumerate(cuts):
            ind = np.where(vals > cut)
            if(np.shape(ind)[1]<1): continue

            nd = str("{:.2f}".format(abs(nds[ic]))).replace('.','_')
            if not h5file:
                outm = outpath+prop+'_cuts_nd'+nd+'_z'+redshift+'.dat'
                tofile = np.column_stack((xgal[ind],ygal[ind],zgal[ind],zgalz[ind]))
                with open(outm,'w') as outf:
                    np.savetxt(outf, tofile, fmt ='%.10e')
                    
            if h5file:
                outm = outpath+prop+'_cuts_nd'+nd+'_z'+redshift+'.h5'
                hf = h5py.File(outm, 'w')

                # Header
                head = hf.create_group('header')
                head.attrs[u'Simulation']    = sim
                head.attrs[u'h0']            = h0
                head.attrs[u'redshift']      = zz
                head.attrs[u'number_density']= nd
                head.attrs[u'units_nd']      = u'log10(nd/(Mpc/h)^-3)'
                head.attrs[u'sel_property']  = prop
                head.attrs[u'units_data']    = u'xgal,ygal,zgal, xgalz, ygalz, zgalz (Mpc/h), log10(mass/Msun/h), log10(sfr/Msun/h/Gyr), log10(lum/h^-2 erg/s), sat(1 sat, 0 cen)'
                print(list(head.attrs.items())) ; sys.exit()

                # Data
                data = hf.create_group('data')
                data.create_dataset('xgal',data=xgal[ind])
                data.create_dataset('ygal',data=ygal[ind])
                data.create_dataset('zgal',data=zgal[ind])
                data.create_dataset('xgalz',data=xgalz[ind])
                data.create_dataset('ygalz',data=ygalz[ind])
                data.create_dataset('zgalz',data=zgalz[ind])
                data.create_dataset('log10mhalo',data=np.log10(mhalo[ind]))
                data.create_dataset('log10mass',data=np.log10(mass[ind]))
                data.create_dataset('log10sfr',data=np.log10(sfr[ind]))
                data.create_dataset('log10lo2_att',data=np.log10(lo2[ind]) + 2*np.log10(h0))
                hf.close()
            print('Output: {}'.format(outm))
    
    f.close()


