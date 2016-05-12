#! /usr/bin/env python
# coding: utf-8

## @file SAGplots_evol.py
## @author Cristian A. Vega Martínez <cnvega(at)fcaglp.unlp.edu.ar>
## @copyright Copyright © 2016 Cristian A. Vega M.
##
## @brief Set of scientific control plots for SAG outputs.
##
## A collection of pre-defined scientific control plots to be
## applied to the SAG hdf5 outputs through the SAGreader module.


import sys
import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pylab as pl
import numpy as np
import math
import h5py

from matplotlib import rcParams

import SAGreader

def SFRvol_z(sagdat, snaplist, outpath, savefile=None, readfile=False, 
             zrange=[0,8], getPlot=False):
   """ Star formation rate density

Routine for generating a star formation rate density evolution in redshift 
plot. An alternative redshift range can be selected.

@param sagdat Input SAGreader.SAGdata or SAGreader.SAGcollection
object data previously loaded
with the corresponding files. Can be replaced by None if 
'readfile' is set. If SAGreader.SAGdata is used, a standard SAG HDF5
output at z=0 including the historical SFR dataset must be loaded.

@param outpath Output folder in which the plot is going to be stored.
The name is chosen automatically by default.

@param savefile (optional) HDF5 file in which the resulting data is 
stored after being calculated. It should not be used together with 
the 'readfile' option.

@param readfile (optional) HDF5 file from which the data is loaded
instead of being read from 'sagdat'. It should not be used together with 
the 'savefile' option.

@param zrange (optional) Redshift range.

@param getPlot (optional) If set to True, the reference to 
matplotlib.pylab is returned by the function at the end instead of
clearing the plot figure.
   """

   print("### Star formation rate density vs redshift")

   datapath = './Data/'
   zmin = min(zrange)
   zmax = max(zrange)

   if not readfile:

      un = sagdat.readUnits()

      if type(sagdat) is SAGreader.SAGdata:
         isReduced = sagdat.reduced
      elif type(sagdat) is SAGreader.SAGcollection:
         isReduced = sagdat.dataList[sagdat.zminidx].reduced

      # 1) Standard outputs:
      if not isReduced:
         # in this case the file should have z=0 galaxies and historical SFR.
         if 0 != sagdat.readAttr('Redshift'):
            print("ERROR: SAG data is not at z=0")
            return
         if 'Histories/SFR' not in sagdat.datasetList():
            print("ERROR: SAG data does not contain 'Histories/SFR' dataset")
            return
            
         BulgeMass = sagdat.readDataset("BulgeMass")
         DiscMass = sagdat.readDataset("DiscMass")

         StellarMass_all = BulgeMass+DiscMass
         del BulgeMass, DiscMass 
         flt_r, flt_c = np.where(StellarMass_all > 0)
         #StellarMass = StellarMass_all[flt_r]
         del StellarMass, flt_c

         sfr_gal_z = sagdat.readDataset("Histories/SFR", idxfilter=flt_r)

         z_all = snaplist.redshifts[:]
         z_flt = np.where((z_min<=z_all)&(z_all<=zmax))[0]

         redshifts = z_all[z_flt]
         sfr = np.zeros(len(z_flt))

         for i, z_idx in enumerate(z_flt):
            sfr[i] = sum(sfr_gal_z[:,z_idx])
         
         sfr *= un.mass.Msun/un.time.yr

      # 2) Reduced outputs (SFR in each snap only, collection needed)
      elif isReduced:
         # loading more than one output could be very time-demanding, so we are
         # not going to filter by stellarmass > 0:
         
         if type(sagdat) is not SAGreader.SAGcollection:
            print("ERROR: A SAGcollection is needed with reduced outputs")
            return

         z_all = np.array(sagdat.redshift[:])
         z_flt = np.where((zmin<=z_all)&(z_all<=zmax))[0]

         redshifts = z_all[z_flt]
         sfr = np.zeros(len(z_flt))

         for i, z_idx in enumerate(z_flt):
            sfr_gal = sagdat.dataList[z_idx].readDataset('SFR')
            sfr[i] = sum(sfr_gal)

         # the reduced files store the SFR in M_sun/h/yr
         sfr /= un.h
  
      sfr /= (sagdat.boxSizeMpc/un.h)**3

   else: # load the SFRd from a file:
      f = h5py.File(readfile, "r")
      sfr = f['SFRdens/sfr'][:] 
      redshifts  = f['SFRdens/z'][:]  
      f.close()

   if savefile:
      f = h5py.File(savefile)
      if 'SFRdens' in f.keys(): del f['SFRdens']
      f['SFRdens/sfr'] = sfr
      f['SFRdens/z']   = redshifts
      f.close()

   # The data for comparison:
   fname = datapath+"SFR_Behroozi_2013.dat"
   ob_x, ob_y, ob_dyh, ob_dyl = np.loadtxt(fname, unpack=True)  

   # Finally, the plot:
   pl.figure()
   pl.errorbar(ob_x, ob_y, yerr=[ob_dyl, ob_dyh], fmt='ok',
               label="Behroozi et al (2013)")
   pl.plot(redshifts, np.log10(sfr), "-r", linewidth=3)
   pl.xlabel(r'$z$', fontsize=16)
   pl.ylabel(r'$\log_{10}({\rm SFR\;density}[{\rm M}_\odot/{\rm yr}/{\rm Mpc}^3])$',
               fontsize=16)
   #pl.xlim((zmin, zmax))
   pl.legend(loc=3, frameon=False, fontsize=12)
   pl.xscale('symlog')
   ticks = [0,0.5,1,1.5,2,4,6,8]
   pl.xticks(ticks, ticks)
   pl.tight_layout()
   if outpath: pl.savefig(outpath+'/SFRdensity_z.eps')
   print("Done!") 
   if getPlot: return pl
   pl.clf()
   return

