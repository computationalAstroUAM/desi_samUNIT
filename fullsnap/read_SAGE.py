# Adam Stevens, 2018
# Functions for reading SAGE data in the format for MultiDark and UNITSIM
import numpy as np
#from pylab import *
import os

"""
These are the properties (and their units) stored inside the file:

SnapNum                     int32       Snapshot number
Type                        int32       Galaxy type (0=central, 1=satellite)
GalaxyIndex                 int64       Unique galaxy identifier for SAGE.  Tracks most massive progenitor
CentralGalaxyIndex          int64       ID for the central galaxy of the same halo
CtreesHaloID                int64       Subhalo index from consistent-trees output.  Note that we have also encoded whether the galaxy is a 'flyby' by setting this ID to be negative if it is a flyby.  The proper consistent-trees output is hence abs(CtreesHaloID).
TreeIndex                   int32       Which (tree)file the galaxy/(sub)halo is stored in
CtreesCentralID             int64       Central halo index from consistent-trees output
mergeType                   int32       0=none; 1=minor merger; 2=major merger; 3=not used; 4=disrupt to intracluster stars
mergeIntoID                 int32       Position in output of galaxy this one merges into (=-1 if not merging)
mergeIntoSnapNum            int32       Snapshot number where the merger occurs
dT                          float32     Time interval of last snapshot [Myr/h]
Pos                       3xfloat32     (x,y,z) coordinates of galaxy [cMpc/h]
Vel                       3xfloat32     (vx,vy,vz) velocity of galaxy [km/s]
Spin                      3xfloat32     Specific angular momentum vector of the (sub)halo [pMpc/h km/s]
Len                         int32       Number of particles in the (sub)halo
Mvir                        float32     M200c for centrals (main haloes) or total mass for subhaloes [10^10 Msun / h]
CentralMvir                 float32     M200c of this galaxy's central halo [10^10 Msun / h]
Rvir                        float32     R200c of the (sub)halo [pMpc/h]
Vvir                        float32     `Virial velocity' of the (sub)halo, i.e. orbital velocity at Rvir [km/s]
Vmax                        float32     Maximum rotation velocity of the (sub)halo [km/s]
VelDisp                     float32     Velocity dispersion of the (sub)halo [km/s]
ColdGas                     float32     Cold (disc) gas mass of the galaxy [10^10 Msun / h]
StellarMass                 float32     Total stellar mass of the galaxy [10^10 Msun / h]
BulgeMass                   float32     Stellar mass of the bulge [10^10 Msun / h]
HotGas                      float32     Hot gas in the (sub)halo [10^10 Msun / h]
EjectedMass                 float32     Ejected hot gas associated with the halo [10^10 Msun / h]
BlackHoleMass               float32     Mass of supermassive black hole [10^10 Msun / h]
IntraClusterStars           float32     Mass of intracluster stars (zero for satellites) [10^10 Msun / h]
MetalsColdGas               float32     Mass of metals in the cold gas [10^10 Msun / h]
MetalsStellarMass           float32     Mass of metals in stars in the galaxy [10^10 Msun / h]
MetalsBulgeMass             float32     Mass of metals in the bulge stars [10^10 Msun / h]
MetalsHotGas                float32     Mass of metals in the hot gas [10^10 Msun / h]
MetalsEjectedMass           float32     Mass of metals in the ejected gas [10^10 Msun / h]
MetalsIntraClusterStars     float32     Mass of metals in the intracluster stars [10^10 Msun / h]
SfrDisk                     float32     Star formation rate of the disk (averaged over the last time-step) [Msun/yr]
SfrBulge                    float32     Star formation rate of the bulge (averaged over the last time-step) [Msun/yr]
SfrDiskZ                    float32     Average metal abundance (ratio of mass in metals to total mass) of gas that formed stars in the disk over the last time-step
SfrBulgeZ                   float32     Average metal abundance (ratio of mass in metals to total mass) of gas that formed stars in the bulge over the last time-step
DiskRadius                  float32     Exponential disk scale radius [pMpc/h]
Cooling                     float32     log10( Net cooling rate [i.e. heating is already subtracted] of hot gas [erg/s] )
Heating                     float32     log10( Heating rate from radio mode AGN that suppresses cooling [erg/s] )
QuasarModeBHaccretionMass   float32     Mass of the black hole the was obtained from quasar mode accretion [10^10 Msun / h] (the remaining mass of the BH comes from radio mode accretion)
TimeOfLastMajorMerger       float32     Look-back time of the last major merger [Myr/h]
TimeOfLastMinorMerger       float32     Look-back time of the last minor merger [Myr/h]
OutflowRate                 float32     Heating rate of cold gas from supernova feedback [Msun / yr]
MeanStarAge                 float32     Mean age (look-back time from z=0) of stars in the galaxy [Myr/h]
infallMvir                  float32     Virial mass of the main halo at the snapshot before this system became a subhalo (only relevant for satellites) [10^10 Msun / h]
infallVvir                  float32     Virial velocity of the main halo at the snapshot before this system became a subhalo [km/s]
infallVmax                  float32     Peak velocity of rotation curve of the main halo at the snapshot before this system became a subhalo [km/s]
"""

def galdtype_multidark():
    Galdesc_full = [
                    ('SnapNum'                      , np.int32),
                    ('Type'                         , np.int32),
                    ('GalaxyIndex'                  , np.int64),
                    ('CentralGalaxyIndex'           , np.int64),
                    ('CtreesHaloID'                 , np.int64),
                    ('TreeIndex'                    , np.int32),
                    ('CtreesCentralID'              , np.int64),
                    ('mergeType'                    , np.int32),
                    ('mergeIntoID'                  , np.int32),
                    ('mergeIntoSnapNum'             , np.int32),
                    ('dT'                           , np.float32),
                    ('Pos'                          , (np.float32, 3)),
                    ('Vel'                          , (np.float32, 3)),
                    ('Spin'                         , (np.float32, 3)),
                    ('Len'                          , np.int32),
                    ('Mvir'                         , np.float32),
                    ('CentralMvir'                  , np.float32),
                    ('Rvir'                         , np.float32),
                    ('Vvir'                         , np.float32),
                    ('Vmax'                         , np.float32),
                    ('VelDisp'                      , np.float32),
                    ('ColdGas'                      , np.float32),
                    ('StellarMass'                  , np.float32),
                    ('BulgeMass'                    , np.float32),
                    ('HotGas'                       , np.float32),
                    ('EjectedMass'                  , np.float32),
                    ('BlackHoleMass'                , np.float32),
                    ('IntraClusterStars'            , np.float32),
                    ('MetalsColdGas'                , np.float32),
                    ('MetalsStellarMass'            , np.float32),
                    ('MetalsBulgeMass'              , np.float32),
                    ('MetalsHotGas'                 , np.float32),
                    ('MetalsEjectedMass'            , np.float32),
                    ('MetalsIntraClusterStars'      , np.float32),
                    ('SfrDisk'                      , np.float32),
                    ('SfrBulge'                     , np.float32),
                    ('SfrDiskZ'                     , np.float32),
                    ('SfrBulgeZ'                    , np.float32),
                    ('DiskRadius'                   , np.float32),
                    ('Cooling'                      , np.float32),
                    ('Heating'                      , np.float32),
                    ('QuasarModeBHaccretionMass'    , np.float32),
                    ('TimeOfLastMajorMerger'        , np.float32),
                    ('TimeOfLastMinorMerger'        , np.float32),
                    ('OutflowRate'                  , np.float32),
                    ('MeanStarAge'                  , np.float32),
                    ('infallMvir'                   , np.float32),
                    ('infallVvir'                   , np.float32),
                    ('infallVmax'                   , np.float32)
                    ]
    names   = [Galdesc_full[i][0] for i in range(len(Galdesc_full))]
    formats = [Galdesc_full[i][1] for i in range(len(Galdesc_full))]
    Galdesc = np.dtype({'names':names, 'formats':formats}, align=True)
    return Galdesc



def sageoutsingle(fname, fields=[]):
    # Read a single SAGE output file, returning all the galaxy data in a record array
    # fname is the full name for the file to read, including its path
    # fields is the list of fields you want to read in.  If empty, will read all fields.
    Galdesc = galdtype_multidark()
    if len(fields)==0:
        fields=list(Galdesc.names)

    fin      = open(fname, 'rb')                                           # Open the file
    Ntrees   = np.fromfile(fin, np.dtype(np.int32),1)                      # Read number of trees in file
    NtotGals = np.fromfile(fin, np.dtype(np.int32),1)[0]                   # Read number of gals in file.
    
    if NtotGals == 0:
        print('no galaxies in file: ',fname)
    
    if (Ntrees != 0):
        GalsPerTree = np.fromfile(fin, np.dtype((np.int32, Ntrees)),1)     # Read the number of gals in each tree
    G = np.fromfile(fin, Galdesc, NtotGals)                                # Read all the galaxy data       
    G = G[fields]                                                          # Reduce to fields of interest
    return G 



def sagesnap(fpre, filelist, fields=[]):
    # Read full SAGE snapshot, going through each file and compiling into 1 array
    # fpre is the name of the file up until the _ before the file number
    # filelist contains all the file numbers you want to read in
    print('  Loading ',fpre,'_',filelist)

    Galdesc = galdtype_multidark()
    if len(fields)==0:
        fields=list(Galdesc.names)

    # create empty storage for all galaxy properties (will be removed again at the end)
#    G = np.empty(1,dtype=Galdesc)
#    G[:][0] = 0.0 
    G = np.zeros(1,dtype=Galdesc)
       
    # extract fields of interest
    G = G[fields]
    
    # loop over file-extensions given in filelist[]
    for i in filelist:
        fname    = fpre+'_'+str(i)
        filesize = os.path.getsize(fname)

        # capture possibly empty files
        if filesize > 0:
            G1 = sageoutsingle(fname, fields)
            G = np.append(G, G1)

    # remove the first (empty) field again
    G = G[:][1:]

    print('    -> found',len(G[:][:]),' galaxies')
    return G


""" Example code for reading in data, feel free to copy:
    fields = ['StellarMass', 'CtreesHaloID', 'Type', 'Mvir'] # specify whichever fields you want to read here (refer to the names in the first function).  Can leave empty
    files = range(100) # list of file numbers to read.  They don't need to be consecutive.
    fpre = '/data2/users/astevens/SAGE_output/UNITSIM1/model_z0.000' # where the data are stored, what the common string in the files is before the last underscore
    G = sagesnap(fpre, files, fields)
    
    # Then do, e.g. G['StellarMass'] to get the stellar masses of the galaxies.
    
    """

