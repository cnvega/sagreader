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

import matplotlib.pyplot as pl
import numpy as np
import math
import h5py

import SAGreader


def set_style(style='book', Hratio=1.0, Wfrac=1.0):
   if style == 'talk':
      size, fsize = 6, 16
   if style == 'book':
      size, fsize = 5.39, 12
   if style == 'mnras':
      size, fsize = 3.32, 8
   if style == 'mnras-fw':
      size, fsize = 6.97, 8

   mpl.rcParams['figure.figsize'] = (size*Wfrac, (size*Wfrac)*Hratio)
   mpl.rcParams['font.family'] = 'serif'
   mpl.rcParams['font.serif'] = ['Times', 'Liberation Serif', 'Times New Roman']
   mpl.rcParams['font.size'] = fsize
   mpl.rcParams['legend.fontsize'] = 'medium'
   mpl.rcParams['legend.frameon'] = False
   mpl.rcParams['text.usetex'] = True

   mpl.rcParams['axes.linewidth'] = 0.5
   try:
      mpl.rcParams['xtick.minor.visible'] = True
      mpl.rcParams['ytick.minor.visible'] = True
   except: pass
   try:
      mpl.rcParams['xtick.top'] = True
      mpl.rcParams['ytick.right'] = True
   except: pass
   try:
      mpl.rcParams['xtick.direction'] = 'in'
      mpl.rcParams['ytick.direction'] = 'in'
   except: pass


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

      if isinstance(sagdat, SAGreader.SAGdata):
         isReduced = sagdat.reduced
         sagdat0 = sagdat
      elif isinstance(sagdat, SAGreader.SAGcollection):
         isReduced = sagdat.dataList[sagdat.zminidx].reduced
         sagdat0 = sagdat.select_redshift(0.0, 0.0)
     
      ### Check if the data is correct:
      if 0 != sagdat0.readAttr('Redshift'):
         print("ERROR: SAG data is not at z=0")
         return
      if 'Histories/SFR' not in sagdat0.datasetList():
         print("ERROR: SAG data does not contain 'Histories/SFR' dataset")
         return

      ### Configuration according to the type of output:
      if isReduced:
         bulgeds, discds = "M_star_bulge", "M_star_disk"
         sfrun = 1.0/un.h    # stored in Msun/yr/h 
      else:
         bulgeds, discds = "BulgeMass", "DiscMass"
         sfrun = un.mass.Msun/un.time.yr  # stored in code units

      z_all = snaplist.redshifts[:]
      z_flt = np.where((zmin<=z_all)&(z_all<=zmax))[0]
      
      redshifts = z_all[z_flt]
      
      sfr = np.zeros(len(z_flt))

      # Using the /Histories/ datasets could demand a lot of RAM, so here
      # we process the loaded files individually:
      for box in range(sagdat0.nfiles):
         if sagdat0.keepOpen: f = sagdat0.dataList[box]
         else: f = h5py.File(sagdat0.filenames[box], "r")
         print("processing file: "+sagdat0.filenames[box])

         sfr_gal_z = f["Histories/SFR"][:]
         
         for i, z_idx in enumerate(z_flt):
            # First, we need to obtain the cosmic SFR:
            sfr[i] += sum(sfr_gal_z[:,z_idx])
            
         if not sagdat0.keepOpen: f.close()
         del sfr_gal_z
      
      # Units and volume: 
      sfr *= sfrun
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
               label="Behroozi et al (2013)", ms=1.5, lw=0.5, zorder=1)
   pl.plot(redshifts, np.log10(sfr), "-r", linewidth=1.6)
   pl.xlabel(r'$z$')
   pl.ylabel(r'$\log (\mathrm{SFR\;density}[\mathrm{M_\odot/yr/Mpc}^3])$')
   pl.xlim((zmin, 6))
   pl.ylim((-3, -0.001))
   pl.tight_layout()
   if outpath: pl.savefig(outpath+'/SFRdensity_z.eps')
   print("Done!") 
   if getPlot: return pl
   pl.clf()
   return


def SFRvol_z_massbin(sagdat, snaplist, outpath, savefile=None, readfile=False, 
                     zrange=[0,8], massbins=[9.7, 10.1, 10.5, 10.9, 11.3], 
                     getPlot=False):
   """ Cosmic star formation rate density for different z=0 stellar mass.

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

   print("### Star formation rate density (stellar mass bins) vs redshift")

   datapath = './Data/'
   zmin = min(zrange)
   zmax = max(zrange)

   massbins = [5]+massbins+[15]    # external bins limits 
   
   if not readfile:

      un = sagdat.readUnits()

      if isinstance(sagdat, SAGreader.SAGdata):
         isReduced = sagdat.reduced
         sagdat0 = sagdat
      elif isinstance(sagdat, SAGreader.SAGcollection):
         isReduced = sagdat.dataList[sagdat.zminidx].reduced
         sagdat0 = sagdat.select_redshift(0.0, 0.0)
     
      ### Check if the data is correct:
      if 0 != sagdat0.readAttr('Redshift'):
         print("ERROR: SAG data is not at z=0")
         return
      if 'Histories/SFR' not in sagdat0.datasetList():
         print("ERROR: SAG data does not contain 'Histories/SFR' dataset")
         return

      ### Configuration according to the type of output:
      if isReduced:
         bulgeds, discds = "M_star_bulge", "M_star_disk"
         sfrun = 1.0/un.h    # stored in Msun/yr/h 
      else:
         bulgeds, discds = "BulgeMass", "DiscMass"
         sfrun = un.mass.Msun/un.time.yr  # stored in code units


      z_all = snaplist.redshifts[:]
      z_flt = np.where((zmin<=z_all)&(z_all<=zmax))[0]
      
      redshifts = z_all[z_flt]
      
      sfr = np.zeros(len(z_flt))
      sfr_bins = np.zeros((len(massbins)-1, len(z_flt)))

      # Using the /Histories/ datasets could demand a lot of RAM, so here
      # we process the loaded files individually:
      for box in range(sagdat0.nfiles):
         if sagdat0.keepOpen: f = sagdat0.dataList[box]
         else: f = h5py.File(sagdat0.filenames[box], "r")
         print("processing file: "+sagdat0.filenames[box])

         BulgeMass, DiscMass = f[bulgeds][:], f[discds][:]

         logsm = np.log10((BulgeMass+DiscMass)*un.mass.Msun/un.h)
         del BulgeMass, DiscMass

         sfr_gal_z = f["Histories/SFR"][:]
         
         for i, z_idx in enumerate(z_flt):
            # First, we need to obtain the cosmic SFR:
            sfr[i] += sum(sfr_gal_z[:,z_idx])
            
            # Now, the contribution for each massbin:
            for b in range(len(massbins)-1):
               mflt = np.where((massbins[b]<=logsm)&(logsm<massbins[b+1]))[0]
               if len(mflt>0):
                  sfr_bins[b,i] += sum(sfr_gal_z[mflt,z_idx])
               del mflt
         
         if not sagdat0.keepOpen: f.close()
         del logsm, sfr_gal_z
      
      # Units and volume: 
      sfr *= sfrun
      sfr_bins *= sfrun
  
      sfr /= (sagdat.boxSizeMpc/un.h)**3
      sfr_bins /= (sagdat.boxSizeMpc/un.h)**3

   else: # load the SFRd from a file:
      f = h5py.File(readfile, "r")
      sfr = f['SFRdens_mb/sfr'][:] 
      redshifts  = f['SFRdens_mb/z'][:]
      sfr_bins   = f['SFRdens_mb/sfr_bins'][:]
      massbins   = f['SFRdens_mb/massbins'][:]
      f.close()

   if savefile:
      f = h5py.File(savefile)
      if 'SFRdens_mb' in f.keys(): del f['SFRdens_mb']
      f['SFRdens_mb/sfr'] = sfr
      f['SFRdens_mb/z']   = redshifts
      f['SFRdens_mb/sfr_bins']   = sfr_bins
      f['SFRdens_mb/massbins']   = massbins
      f.close()

   # The data for comparison:
   fname = datapath+"SFR_Behroozi_2013.dat"
   ob_x, ob_y, ob_dyh, ob_dyl = np.loadtxt(fname, unpack=True)  
   
   ### Load another model:
   #ff = h5py.File("plots/SAG-7.130_aF1.3/data.h5", "r")
   #xx = ff['SFRdens/z'][:]
   #yy = ff['SFRdens/sfr'][:]
   #ff.close()

   # Finally, the plot:
   pl.figure()
   pl.errorbar(ob_x, ob_y, yerr=[ob_dyl, ob_dyh], fmt='ok',
               label="Behroozi et al (2013)", ms=1.5, lw=0.5, zorder=1)
   #pl.plot(xx, np.log10(yy), "--b", linewidth=1.6, 
   #        label=r"SAG, $\beta_\mathrm{F}=1.3$")
   pl.plot(redshifts, np.log10(sfr), "-r", linewidth=1.6, 
           label="SAG")
   lshape = ['-','--','-.',(0,(4,2,1,2,1,2)),(0,(4,1,1,1,1,1,1,1)),':']
   lshape.reverse()
   for i in range(len(massbins)-1):
      if i==0:
         label = r"$ \log M_\star < "+str(massbins[1])+"$"
      elif i==len(massbins)-2:
         label = r"$"+str(massbins[-2])+" \leq \log M_\star $"
      else:
         label = r"$"+str(massbins[i])+" \leq \log M_\star < "+str(massbins[i+1])+"$"
      
      pl.plot(redshifts, np.log10(sfr_bins[i]), "k", linestyle=lshape[i], lw=1,
              label=label)
   pl.xlabel(r'$z$')
   pl.ylabel(r'$\log (\mathrm{SFR\;density}[\mathrm{M_\odot/yr/Mpc}^3])$')
   pl.xlim((zmin, 6))
   pl.ylim((-3, -0.001))
   #pl.legend(loc=3, frameon=False)
   leg = pl.legend(fontsize='x-small', frameon=True, edgecolor='k', 
         handlelength=2.5, labelspacing=0.2)
   leg.get_frame().set_linewidth(0.5)
   
   pl.tight_layout()
   if outpath: pl.savefig(outpath+'/SFRdensity_z_mb.eps')
   print("Done!") 
   if getPlot: return pl
   pl.clf()
   return

