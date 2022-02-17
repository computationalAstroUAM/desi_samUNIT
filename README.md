To-do:
[] Measure bias in z-space and compare

Maybe:
[] hdf5 file with further properties
[] HODs of the samples

In this folder there are codes to generate a [DESI-like catalogues](https://desi.lbl.gov/trac/wiki/Clustering/MockChallenge/make_galaxy) from the SAGE+UNIT catalogue.

## Programs

**lo2f.py** To plot the model [OII] luminosity function at z~1, compared with observations.

**cumu.py** To plot and create a file with the cumulative function on the indicated property and find the cuts on it that will give the target DESI number density.

**cumu_multi.py** Cumulative function for only those galaxies not selected as 'elgs' with the cuts done above.

**find_cuts.py** Find the cuts in a given property to get a target number density

**create_sel.py** Use the found cuts to create a subsample.

**z_zrsd.py** Plot zrsd vs z to check the created subsamples.

**pklin.py** Calculate the matter power spectrum P(k,z) from CAMB (with and without non linear corrections)

**pkg.py** Calculate P(k) for the sample galaxies using nbodykit

**bias_r_pk.py** Obtain the bias in real space for the different selections.

**bias_z_pk.py** Obtain the bias in z-space for the different selections. WORK IN PROGRESS

## Subfolders

**fullsnap** Produce a file to be shared, from a full snapshot of SAGE

## Specifications for DESI samples

Specifications following the mocks tabulated [here](https://desi.lbl.gov/trac/wiki/Clustering/MockChallenge/post-recon-BAO/stage2).

###UNITSIM
From DESI (https://desi.lbl.gov/trac/wiki/Clustering/MockChallenge/post-recon-BAO/stage2):

*ELGs* mocks at redshift = 0.9873 (snap 97), densities={25e-4,20e-4,5.5e-4}(Mpc/h)^-3.

*LRGs* mocks at z=0.9873, density=4.4e-4 (Mpc/h)^-3 and z=0.7400, density=5.5e-4 (Mpc/h)^-3.

In taurus, emission lines for the HighRes UNITSIM in [data6](/data6/users/aknebe/Projects/UNITSIM/ELGs_DESI/):

* Flux [erg/s/cm^2]
* Luminosity [erg/s]

## Output

### Format
The output should be a text file with (following [DESI wiki](https://desi.lbl.gov/trac/wiki/CosmoSimsWG/DESI_mocks)):

{x, y, z, z_RSD} in units of Mpc/h

z_RSD is including the position correction due to redshift distortion in the z-axis. 

z_RSD = z + vz*(1+redshift)/H,

where H=100*sqrt(Omega_m*(1+redshift)3 + Omega_Lambda)


Using awk from default ASCII output files:
```
awk '{print $1 " " $2 " " $3 " " $6}' inputfile > outfile
```

### Location at NERSC

Once the catalogues are tested, they can go to NERSC:
```
/global/project/projectdirs/desi/mocks/UNIT/SAGE_[ELG, LRG, QSO]/z*
```
