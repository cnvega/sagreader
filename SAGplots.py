#! /usr/bin/env python
# coding: utf-8

## @file SAGplots.py
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
   mpl.rcParams['font.serif'] = ['Liberation Serif', 'Times New Roman']
   mpl.rcParams['font.size'] = fsize
   mpl.rcParams['legend.fontsize'] = 'medium'
   mpl.rcParams['legend.frameon'] = False
   mpl.rcParams['xtick.minor.visible'] = True
   mpl.rcParams['ytick.minor.visible'] = True



def SMF(sagdat, outpath, savefile=None, readfile=False, redshift=0,
        getPlot=False):
   """ Stellar Mass function

Routine for generating a stellar mass function plot. It can be used
indistinctly for the redshifts 0, 1, 2 and 3.

@param sagdat Input SAGreader.SAGdata object data previously loaded
with the corresponding redshift. Can be replaced by None if 
'readfile' is set.

@param outpath Output folder in which the plot is going to be stored.
The name is chosen automatically by default.

@param savefile (optional) HDF5 file in which the resulting data is 
stored after being calculated. It should not be used together with 
the 'readfile' option.

@param readfile (optional) HDF5 file from which the data is loaded
instead of being read from 'sagdat'. It should not be used together with 
the 'savefile' option.

@param redshift (optional) Integer value indicating the redshift of 
the loaded data. Is can be either 0, 1, 2 or 3, and is used for loading
the corresponding observational data.

@param getPlot (optional) If set to True, the reference to 
matplotlib.pylab is returned by the function at the end instead of
clearing the plot figure.
   """
   print("### Stellar Mass Function")

   datapath = './Data/'
  
   if redshift in [0,1,2,3]:
      z = redshift
   else:
      print("Redshift not available!")
      return 
   fname = datapath+"smfz"+str(z)+"_Henriques.dat"
   x1, x2, ob_y, ob_e = np.loadtxt(fname, unpack=True)
   ob_x = (x1 + x2)/2.
   ob_eu = np.log10((ob_y+ob_e)/ob_y)
   ob_ed = np.log10(ob_y/(ob_y-ob_e))

   if not readfile:

      if sagdat.reduced:
         BulgeMass = sagdat.readDataset("M_star_bulge")
         DiscMass = sagdat.readDataset("M_star_disk")

      else:
         BulgeMass = sagdat.readDataset("BulgeMass")
         DiscMass = sagdat.readDataset("DiscMass")

      # Units:
      un = sagdat.readUnits()
      BulgeMass *= un.mass.Msun
      DiscMass *= un.mass.Msun

      StellarMass_all = BulgeMass+DiscMass
      del BulgeMass, DiscMass 
      StellarMass = StellarMass_all[StellarMass_all > 0]
      del StellarMass_all

      # Correct by h:
      StellarMass /= un.h
      # log scale: 
      logSM = np.log10(StellarMass) 
      del StellarMass

      Vol = sagdat.boxSizeMpc**3

      # And the counting:
      phi_i, bins = np.histogram(logSM, bins=45, range=[6,14])

      phi = phi_i.astype(float)/Vol/(bins[1]-bins[0])

   else: # load the histograms from a file:
      f = h5py.File(readfile, "r")
      bins = f['SMF/z'+str(z)+'/bins'][:] 
      phi  = f['SMF/z'+str(z)+'/phi'][:]  
      f.close()

   if savefile:
      f = h5py.File(savefile)
      if 'SMF/z'+str(z) in f.keys(): del f['SMF/z'+str(z)]
      f['SMF/z'+str(z)+'/bins'] = bins
      f['SMF/z'+str(z)+'/phi']   = phi
      f.close()

   # Finally, the plot:
   pl.figure()
   x = (bins[1:]+bins[:-1])/2.0
   pl.plot(x[phi>0], np.log10(phi[phi>0]), '-k', lw=2, label='SAG')
   pl.errorbar(ob_x, np.log10(ob_y) , yerr=[ob_ed, ob_eu], fmt='^b',
              label="Henriques et al (2013), $z="+str(z)+"$")
   pl.xlabel(r'$\log_{10}(M_\star[{\rm M}_\odot])$')
   pl.ylabel(r'$\log_{10}(\Phi[{\rm h}^{3} {\rm Mpc}^{-3} /\log_{10} M_\star])$')
   pl.xlim((7,13))
   pl.legend(loc=3, frameon=False)
   #pl.axis('scaled')
   pl.tight_layout()
   pl.savefig(outpath+'/SMF_z'+str(z)+'.eps')
   print("Done!")
   if getPlot: return pl
   pl.clf()
   return 


def FracMorph(sagdat, outpath, savefile=None, readfile=False,
              threshSE=0.85, threshIS=0.0, SEDmagfilter='NONE',
              getPlot=False):
   """ Galaxy morphological fraction 

Routine for generating galaxy morphological fraction plot at z=0,
by comparing the mass fractions of galaxy components and applying 
threshold to classify them. A 'V' magnitude filter is tried first
for selecting the galaxies before the calculation, by loading two known 
datasets or the defined by the SEDmagfilter option. If none of them is 
found, all the galaxies are included in the selection. 

@param sagdat Input SAGreader.SAGdata or SAGreader.SAGcollection
object data previously loaded
with the corresponding files. Can be replaced by None if 
'readfile' is set.

@param outpath Output folder in which the plot is going to be stored.
The name is chosen automatically by default.

@param savefile (optional) HDF5 file in which the resulting data is 
stored after being calculated. It should not be used together with 
the 'readfile' option.

@param readfile (optional) HDF5 file from which the data is loaded
instead of being read from 'sagdat'. It should not be used together with 
the 'savefile' option.

@param threshSE (optional) Bulge to total mass ratio threshold to
classify elliptical and disc galaxies. The default value (0.85) is used
if omitted.

@param threshIS (optional) Bulge to total mass ratio threshold to
classify irregular and disc galaxies. The default value (0.0) is used
if omitted. 

@param SEDmagfilter (optional) Alternative magnitude filter (inside the
'SED/Magnitudes' group) for being used for filtering the galaxies.

@param getPlot (optional) If set to True, the reference to 
matplotlib.pylab is returned by the function at the end instead of
clearing the plot figure.
   """ 
   print("### Morphological Fraction")

   nbins = 40
   mrange = [7,13]
   
   datapath = './Data/'

   if not readfile:
      
      # First we check if there are an available Magnitude for pre-selecting the 
      # galaxies:

      if 'Magnitudes/Mag_V_dust1' in sagdat.datasetList():
         MagV = sagdat.readDataset('Magnitudes/Mag_V_dust1')
         mfilt = MagV < 0
         del MagV
      if 'Magnitudes/Mag_V_dust' in sagdat.datasetList():
         MagV = sagdat.readDataset('Magnitudes/Mag_V_dust')
         mfilt = MagV < 0
         del MagV
      elif 'SED/Magnitudes/Mag_ext_id14_tot_r' in sagdat.datasetList():
         MagV = sagdat.readDataset('SED/Magnitudes/Mag_ext_id14_tot_r')
         mfilt = MagV < 0
         del MagV
      elif 'SED/Magnitudes/'+SEDmagfilter in sagdat.datasetList():
         MagV = sagdat.readDataset('SED/Magnitudes/'+SEDmagfilter)
         mfilt = MagV < 0
         del MagV
      else:
         print("Warning: No magnitude filter applied in the selection.")
         mfilt = []

      if sagdat.reduced:
         BulgeMass = sagdat.readDataset("M_star_bulge", idxfilter=mfilt)
         DiscMass = sagdat.readDataset("M_star_disk", idxfilter=mfilt)

      else:
         BulgeMass = sagdat.readDataset("BulgeMass", idxfilter=mfilt)
         DiscMass = sagdat.readDataset("DiscMass", idxfilter=mfilt)

      # Units:
      un = sagdat.readUnits()
      BulgeMass *= un.mass.Msun
      DiscMass *= un.mass.Msun
      Hubble_h = un.h

      wmass = (BulgeMass + DiscMass) > 0

      logSM = np.log10(BulgeMass[wmass] + DiscMass[wmass])
      ratio = BulgeMass[wmass]/(DiscMass[wmass]+BulgeMass[wmass])
      del BulgeMass, DiscMass

      # Now the counting:
      tot, bin_tot = np.histogram(logSM, bins=nbins, range=mrange)
      irr, bin_irr = np.histogram(logSM[ratio<=threshIS], bins=nbins, range=mrange)
      sp, bin_sp  = np.histogram(logSM[(threshIS<ratio)&(ratio<=threshSE)],
                                 bins=nbins, range=mrange)
      el, bin_el  = np.histogram(logSM[ratio>threshSE], bins=nbins, range=mrange)

      del logSM, ratio

   else: # read histograms from file:
      f = h5py.File(readfile, "r")
      tot     = f['Morph/total_counts'][:]
      bin_tot = f['Morph/bins'][:]
      irr = f['Morph/irr'][:]
      sp  = f['Morph/sp'][:]
      el  = f['Morph/el'][:]
      Hubble_h = f['Morph/h'].value
      f.close()

   if savefile:
      f = h5py.File(savefile)
      if 'Morph' in f.keys(): del f['Morph']
      f['Morph/total_counts'] = tot
      f['Morph/bins'] = bin_tot
      f['Morph/irr'] = irr 
      f['Morph/sp']  = sp  
      f['Morph/el']  = el  
      f['Morph/h']  = Hubble_h
      f.close()

   x = ((bin_tot[:-1]+bin_tot[1:])/2.)[tot > 0]
   f_sp  = sp.astype(np.float)[tot>0] / tot.astype(np.float)[tot>0]
   f_el  = el.astype(np.float)[tot>0] / tot.astype(np.float)[tot>0]
   f_irr = irr.astype(np.float)[tot>0] / tot.astype(np.float)[tot>0]

   # Loading the observations:
   obs_sp = np.loadtxt(datapath+"conselice2006_Sp.dat", unpack=True)
   obs_el = np.loadtxt(datapath+"conselice2006_Ell.dat", unpack=True)
   obs_irr = np.loadtxt(datapath+"conselice2006_Irr.dat", unpack=True)

   # Finally, the plot:
   pl.figure()
   pl.plot(x, f_sp, "b--", label="Spirals", linewidth=2)
   pl.plot(x, f_el, "r-.", label="Ellipticals", linewidth=2)
   pl.plot(x, f_irr, "k:", label="Irregulars", linewidth=2)
   pl.legend(frameon=False, loc='upper left')

   pl.errorbar(obs_sp[0]+np.log10(Hubble_h), obs_sp[1], yerr=obs_sp[2], fmt="db")
   pl.errorbar(obs_el[0]+np.log10(Hubble_h), obs_el[1], yerr=obs_el[2], fmt="or")
   pl.errorbar(obs_irr[0]+np.log10(Hubble_h), obs_irr[1], yerr=obs_irr[2], fmt="^k")

   pl.xlabel(r"$\log_{10}(M_\star [h^{-1} {\rm M}_\odot])$")
   pl.ylabel("Fraction")
   pl.xlim((8.,12.))
   pl.xticks(np.arange(8,13))
   pl.ylim((-0.05, 1.05))
   pl.tight_layout()
   pl.savefig(outpath+'/FracMorph.eps')
   print("Done!")
   if getPlot: return pl
   pl.clf()
   

def BHBulge(sagdat, outpath, savefile=None, readfile=False, 
            SEDmagfilter='NONE', getPlot=False):
   """ Black hole mass vs bulge mass relationship.  

Routine for generating the super-massive black hole mass versus
bulge mass relationship plot at z=0.
A 'V' magnitude filter is tried first
for selecting the galaxies before the calculation, by loading two known 
datasets or the defined by the SEDmagfilter option. If none of them is 
found, all the galaxies are included in the selection. 

@param sagdat Input SAGreader.SAGdata or SAGreader.SAGcollection
object data previously loaded
with the corresponding files. Can be replaced by None if 
'readfile' is set.

@param outpath Output folder in which the plot is going to be stored.
The name is chosen automatically by default.

@param savefile (optional) HDF5 file in which the resulting data is 
stored after being calculated. It should not be used together with 
the 'readfile' option.

@param readfile (optional) HDF5 file from which the data is loaded
instead of being read from 'sagdat'. It should not be used together with 
the 'savefile' option.

@param SEDmagfilter (optional) Alternative magnitude filter (inside the
'SED/Magnitudes' group) for being used for filtering the galaxies.

@param getPlot (optional) If set to True, the reference to 
matplotlib.pylab is returned by the function at the end instead of
clearing the plot figure.
   """ 
   print("### Black hole vs Bulge mass relationship")

   from matplotlib import colors, ticker, cm
   from numpy import ma

   datapath = './Data/'

   xr = [8,13]
   yr = [5,11]

   if not readfile:
      if 'Magnitudes/Mag_V_dust1' in sagdat.datasetList():
         MagV = sagdat.readDataset('Magnitudes/Mag_V_dust1')
         mfilt = MagV < 0
         del MagV
      if 'Magnitudes/Mag_V_dust' in sagdat.datasetList():
         MagV = sagdat.readDataset('Magnitudes/Mag_V_dust')
         mfilt = MagV < 0
         del MagV
      elif 'SED/Magnitudes/Mag_ext_id14_tot_r' in sagdat.datasetList():
         MagV = sagdat.readDataset('SED/Magnitudes/Mag_ext_id14_tot_r')
         mfilt = MagV < 0
         del MagV
      elif 'SED/Magnitudes/'+SEDmagfilter in sagdat.datasetList():
         MagV = sagdat.readDataset('SED/Magnitudes/'+SEDmagfilter)
         mfilt = MagV < 0
         del MagV
      else:
         print("Warning: No magnitude filter applied in the selection.")
         mfilt = []


      if sagdat.reduced:
         BulgeMass = sagdat.readDataset("M_star_bulge", idxfilter=mfilt)
         BHMass = sagdat.readDataset("Mbh", idxfilter=mfilt)

      else:
         BulgeMass = sagdat.readDataset("BulgeMass", idxfilter=mfilt)
         BHMass = sagdat.readDataset("AGN/BHMass", idxfilter=mfilt)

      # Units:
      un = sagdat.readUnits()
      BulgeMass *= un.mass.Msun/un.h
      BHMass *= un.mass.Msun/un.h
      
      nz = (BulgeMass>0)&(BHMass>0)

      H, xedges, yedges = np.histogram2d(np.log10(BulgeMass[nz]),
                                np.log10(BHMass[nz]), range=[xr,yr], bins=80)
      del BulgeMass, BHMass

   else: #read histograms from file:
      f = h5py.File(readfile, "r")
      H = f['BHBulge/histogram'][:]
      xedges = f['BHBulge/Xedges_bulge'][:]
      yedges = f['BHBulge/Yedges_bh'][:]
   
   if savefile:
      f = h5py.File(savefile)
      if 'BHBulge' in f.keys(): del f['BHBulge']
      f['BHBulge/histogram']    = H
      f['BHBulge/Xedges_bulge'] = xedges
      f['BHBulge/Yedges_bh']    = yedges


   # to avoid the warning with 0 in log.
   H = ma.masked_where(H<=0, H)
   #x = (xedges[1:]+xedges[:-1])/2
   #y = (yedges[1:]+yedges[:-1])/2

   # Observations:
   obK = np.loadtxt(datapath+"BHB_Kormendy2013.dat", unpack=True)
   obM = np.loadtxt(datapath+"BHB_McConnell2013.dat", unpack=True)

   # and the plot:
   fig, ax = pl.subplots(1,1)
   
   ax.errorbar(obK[0], obK[2], xerr=obK[1], yerr=[obK[4], obK[3]], 
               fmt="ob", label="Kormendy & Ho (2013)", ms=4)
   ax.errorbar(obM[0], obM[2], xerr=obM[1], yerr=[obM[4], obM[3]], 
               fmt="^g", label="McConnell & Ma (2013)", ms=4)

   # This is just for not including the error bars in the legend
   handles, labels = ax.get_legend_handles_labels()
   handles = [h[0] for h in handles]
   ax.legend(handles, labels, frameon=False, loc="upper left",
             numpoints=1)

   #cs = ax.contourf(x, y, H, vmin=1, locator=ticker.LogLocator(), cmap=cm.Reds)
   cs = ax.imshow(H.T, origin="low", cmap=cm.YlOrBr, norm=colors.LogNorm(vmin=1),
                  extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]], 
                  aspect='auto', interpolation='nearest')
   fig.colorbar(cs)


   ax.set_xlabel(r"$\log_{10}(M_{\rm Bulge} [{\rm M}_\odot])$")
   ax.set_ylabel(r"$\log_{10}(M_{\rm BH} [{\rm M}_\odot])$")
   
   ax.set_xlim(xr)
   ax.set_ylim(yr)

   pl.tight_layout()

   fig.savefig(outpath+"/BHBulge.eps")
   print("Done!")
   if getPlot: return pl
   pl.clf()
   

def TullyFisher(sagdat, outpath, savefile=None, readfile=False, 
                SEDmagfilter='NONE', getPlot=False):
   """ Tully-Fisher relationship.  

Routine for generating the Tully-Fisher relationship plot at z=0.
A 'V' magnitude filter is tried first
for selecting the galaxies before the calculation, by loading two known 
datasets or the defined by the SEDmagfilter option. If none of them is 
found, all the galaxies are included in the selection. 

@param sagdat Input SAGreader.SAGdata or SAGreader.SAGcollection
object data previously loaded
with the corresponding files. Can be replaced by None if 
'readfile' is set.

@param outpath Output folder in which the plot is going to be stored.
The name is chosen automatically by default.

@param savefile (optional) HDF5 file in which the resulting data is 
stored after being calculated. It should not be used together with 
the 'readfile' option.

@param readfile (optional) HDF5 file from which the data is loaded
instead of being read from 'sagdat'. It should not be used together with 
the 'savefile' option.

@param SEDmagfilter (optional) Alternative magnitude filter (inside the
'SED/Magnitudes' group) for being used for filtering the galaxies.

@param getPlot (optional) If set to True, the reference to 
matplotlib.pylab is returned by the function at the end instead of
clearing the plot figure.
   """ 
   print("### Tully Fisher relationship")
   from matplotlib import colors, ticker, cm
   from numpy import ma

   datapath = './Data/'

   xr = [1.35,   2.8]
   yr = [-23.5, -11]

   if not readfile:
      if sagdat.reduced:
         print("ERROR: VcFlat dataset not included in reduced HDF5 outputs")
         return
         
      # Filters: V mag and galaxy type:
 
      galType = sagdat.readDataset("Type")

      if 'Magnitudes/Mag_V_dust1' in sagdat.datasetList():
         MagV = sagdat.readDataset('Magnitudes/Mag_V_dust1')
         flt = (MagV < 0)&(0 == galType)
         del MagV, galType
      if 'Magnitudes/Mag_V_dust' in sagdat.datasetList():
         MagV = sagdat.readDataset('Magnitudes/Mag_V_dust')
         flt = (MagV < 0)&(0 == galType)
         del MagV, galType
      elif 'SED/Magnitudes/Mag_ext_id14_tot_r' in sagdat.datasetList():
         MagV = sagdat.readDataset('SED/Magnitudes/Mag_ext_id14_tot_r')
         flt = (MagV < 0)&(0 == galType)
         del MagV, galType
      elif 'SED/Magnitudes/'+SEDmagfilter in sagdat.datasetList():
         MagV = sagdat.readDataset('SED/Magnitudes/'+SEDmagfilter)
         flt = (MagV < 0)&(0 == galType)
         del MagV, galType
      else:
         print("Warning: No magnitude V filter applied in the selection.")
         flt = (0 == galType)
      
      # Mag_rS
      if 'Magnitudes/Bulge/Mag_rSb' in sagdat.datasetList():
         if 'Magnitudes/Mag_rS_dust1' in sagdat.datasetList():
            Mag_rS = sagdat.readDataset('Magnitudes/Mag_rS_dust1', idxfilter=flt)
         elif 'Magnitudes/Mag_rS_dust' in sagdat.datasetList():
            Mag_rS = sagdat.readDataset('Magnitudes/Mag_rS_dust', idxfilter=flt)
         else:
            print("ERROR: SDSS-r total magnitudes not found in output")
            return
         
         Mag_rSb = sagdat.readDataset('Magnitudes/Bulge/Mag_rSb', idxfilter=flt)
         
      elif 'SED/Magnitudes/Mag_ext_id121_bulge_r' in sagdat.datasetList():
         Mag_rS = sagdat.readDataset('SED/Magnitudes/Mag_ext_id121_tot_r', idxfilter=flt)
         Mag_rSb = sagdat.readDataset('SED/Magnitudes/Mag_ext_id121_bulge_r', idxfilter=flt)
      else:
         print("ERROR: SDSS-r bulge magnitudes not found in output")
         return

      VcFlat = sagdat.readDataset('VcFlat', idxfilter=flt)

      # Units:
      un = sagdat.readUnits()

      # Apply the same filter as Guo et al. (2011), selecting the galaxies for which
      # the r-band magnitude of the bulge is at least 1.5 mag fainter than that of 
      # the galaxy as a whole.
      deltaMag = 1.5
      guoflt = (np.abs(Mag_rS - Mag_rSb) > deltaMag)
      del Mag_rSb

      Magr = Mag_rS[guoflt] - 5.*np.log10(un.h)
      del Mag_rS

      Vrot = VcFlat[guoflt]*un.velocity.km_per_s
      del VcFlat

      H, xedges, yedges = np.histogram2d(np.log10(Vrot), Magr, range=[xr,yr], bins=80)

   else: #read histograms from file:
      f = h5py.File(readfile, "r")
      H = f['TullyFisher/histogram']
      xedges = f['TullyFisher/Xedges_Vrot']
      yedges = f['TullyFisher/Yedges_Magr']
   
   if savefile:
      f = h5py.File(savefile)
      if 'TullyFisher' in f.keys(): del f['TullyFisher']
      f['TullyFisher/histogram'] = H
      f['TullyFisher/Xedges_Vrot'] = xedges
      f['TullyFisher/Yedges_Magr'] = yedges


   # to avoid the warning with 0 in log.
   H = ma.masked_where(H<=0, H)
   #x = (xedges[1:]+xedges[:-1])/2
   #y = (yedges[1:]+yedges[:-1])/2
   
   # Observations:
   obs = np.loadtxt(datapath+'TF_rband_guo2011.dat', unpack=True)
   
   # and the plot:
   fig, ax = pl.subplots(1,1)
   
   ax.errorbar(obs[0], obs[2], xerr=obs[1], yerr=obs[3], 
               fmt="om", label="Blanton et al. (2008), Springob et al (2007)")

   # This is just for not including the error bars in the legend
   handles, labels = ax.get_legend_handles_labels()
   handles = [h[0] for h in handles]
   ax.legend(handles, labels, frameon=False, loc="upper left",
             numpoints=1)

   #cs = ax.contourf(x, y, H, vmin=1, locator=ticker.LogLocator(), cmap=cm.Greens)
   cs = ax.imshow(H.T, origin='lower', 
                  cmap=cm.Greens, norm=colors.LogNorm(vmin=1),
                  extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]], 
                  aspect='auto', interpolation='nearest')
   fig.colorbar(cs)

   #ax.plot(np.log10(Vrot), Magr, ".k", markersize=0.1)

   ax.set_xlabel(r"$\log_{10}(V_{\rm flat} [{\rm km/s}])$")
   ax.set_ylabel(r"$M_{\rm r} - 5 \log_{10} h$")
   
   ax.set_xlim(xr)
   ax.set_ylim([yr[1],yr[0]])

   pl.tight_layout()

   fig.savefig(outpath+"/TullyFisher.eps")
   print("Done!")
   if getPlot: return pl
   pl.clf()
   
   
def CMD(sagdat, outpath, savefile=None, readfile=False, 
        masscut=0., divpar=(1.8, 18.7), getPlot=False):
   """ Color-magnitud diagram.  

Routine for generating the color-magnitud distribution plot at z=0.

@param sagdat Input SAGreader.SAGdata or SAGreader.SAGcollection
object data previously loaded
with the corresponding files. Can be replaced by None if 
'readfile' is set.

@param outpath Output folder in which the plot is going to be stored.
The name is chosen automatically by default.

@param savefile (optional) HDF5 file in which the resulting data is 
stored after being calculated. It should not be used together with 
the 'readfile' option.

@param readfile (optional) HDF5 file from which the data is loaded
instead of being read from 'sagdat'. It should not be used together with 
the 'savefile' option.

@param masscut (optional) Mass minimum limit in Msun/h.

@param divpar (optional) Free parameters for the division between 
the red and blue galaxy population.

@param getPlot (optional) If set to True, the reference to 
matplotlib.pylab is returned by the function at the end instead of
clearing the plot figure.
   """ 
   print("### Color-magnitud diagram")
   from matplotlib import colors, ticker, cm
   from numpy import ma

   datapath = './Data/'

   xr = [-23,  -15]
   yr = [0.5, 3.3]

   if not readfile:
      
      # let's apply a mass filter
      if sagdat.reduced:
         BulgeMass = sagdat.readDataset("M_star_bulge")
         DiscMass = sagdat.readDataset("M_star_disk")

      else:
         BulgeMass = sagdat.readDataset("BulgeMass")
         DiscMass = sagdat.readDataset("DiscMass")

      # Units:
      un = sagdat.readUnits()
      
      StellarMass = (BulgeMass+DiscMass)*un.mass.Msun
      del BulgeMass, DiscMass 
      massflt = StellarMass > masscut
      del StellarMass
 
      if 'Magnitudes/Mag_uS_dust1' in sagdat.datasetList():
         uMag = sagdat.readDataset('Magnitudes/Mag_uS_dust1', idxfilter=massflt)
         rMag = sagdat.readDataset('Magnitudes/Mag_rS_dust1', idxfilter=massflt)
      elif 'Magnitudes/Mag_uS_dust' in sagdat.datasetList():
         uMag = sagdat.readDataset('Magnitudes/Mag_uS_dust', idxfilter=massflt)
         rMag = sagdat.readDataset('Magnitudes/Mag_rS_dust', idxfilter=massflt)
      elif 'SED/Magnitudes/Mag_ext_id119_tot_r' in sagdat.datasetList():
         uMag = sagdat.readDataset('SED/Magnitudes/Mag_ext_id119_tot_r')
         rMag = sagdat.readDataset('SED/Magnitudes/Mag_ext_id121_tot_r')
      else:
         print("ERROR: SDSS-u,r total magnitudes not found in output")
         return
      
      xplot = rMag - 5*np.log10(un.h/0.7)
      yplot = uMag - rMag

      H, xedges, yedges = np.histogram2d(xplot, yplot, range=[xr,yr], bins=80)

   else: #read histograms from file:
      f = h5py.File(readfile, "r")
      H = f['CMD/histogram']
      xedges = f['CMD/Xedges']
      yedges = f['CMD/Yedges']
   
   if savefile:
      f = h5py.File(savefile)
      if 'CMD' in f.keys(): del f['CMD']
      f['CMD/histogram'] = H
      f['CMD/Xedges'] = xedges
      f['CMD/Yedges'] = yedges


   # to avoid the warning with 0 in log.
   H = ma.masked_where(H<=0, H)
   
   # and the plot:
   fig, ax = pl.subplots(1,1)
   
   cs = ax.imshow(H.T, origin='lower', 
                  cmap=cm.YlGn, norm=colors.LogNorm(vmin=1),
                  extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]], 
                  aspect='auto', interpolation='nearest')
   fig.colorbar(cs)

   xdiv = np.arange(xr[0], xr[1], 0.01)
   ydiv = divpar[0] - 0.2444*np.tanh((xdiv+divpar[1])/1.09)

   ax.plot(xdiv, ydiv, 'm-', lw=1)

   ax.set_xticks([-22, -20, -18, -16])
   ax.set_xlabel(r"$r - 5\log_{10}(h_{70})$")
   ax.set_ylabel(r"$u - r$")
   
   ax.set_xlim(xr)
   ax.set_ylim(yr)

   pl.tight_layout()

   fig.savefig(outpath+"/CMD_uur.eps")
   print("Done!")
   if getPlot: return pl
   pl.clf()
  

def RedFraction(sagdat, outpath, savefile=None, readfile=False, 
        masscut=0., divpar=(1.8, 18.7), getPlot=False):
   """ Passive/red fraction.  

Routine for generating the fraction of red galaxies.

@param sagdat Input SAGreader.SAGdata or SAGreader.SAGcollection
object data previously loaded
with the corresponding files. Can be replaced by None if 
'readfile' is set.

@param outpath Output folder in which the plot is going to be stored.
The name is chosen automatically by default.

@param savefile (optional) HDF5 file in which the resulting data is 
stored after being calculated. It should not be used together with 
the 'readfile' option.

@param readfile (optional) HDF5 file from which the data is loaded
instead of being read from 'sagdat'. It should not be used together with 
the 'savefile' option.

@param masscut (optional) Mass minimum limit in Msun/h.

@param divpar (optional) Free parameters for the division between 
the red and blue galaxy population.

@param getPlot (optional) If set to True, the reference to 
matplotlib.pylab is returned by the function at the end instead of
clearing the plot figure.
   """ 
   print("### Red/passive fraction")
   from matplotlib import colors, ticker, cm
   from numpy import ma

   datapath = './Data/'

   #xr = [-23,  -15]
   #yr = [0.5, 3.3]

   if not readfile:
      
      # loading the stellar mass
      if sagdat.reduced:
         BulgeMass = sagdat.readDataset("M_star_bulge")
         DiscMass = sagdat.readDataset("M_star_disk")

      else:
         BulgeMass = sagdat.readDataset("BulgeMass")
         DiscMass = sagdat.readDataset("DiscMass")

      # Units:
      un = sagdat.readUnits()
      
      StellarMass_all = (BulgeMass+DiscMass)*un.mass.Msun
      del BulgeMass, DiscMass 
      massflt = StellarMass_all > masscut
      StellarMass = StellarMass_all[massflt]
      del StellarMass_all
 
      if 'Magnitudes/Mag_uS_dust1' in sagdat.datasetList():
         uMag = sagdat.readDataset('Magnitudes/Mag_uS_dust1', idxfilter=massflt)
         rMag = sagdat.readDataset('Magnitudes/Mag_rS_dust1', idxfilter=massflt)
      elif 'Magnitudes/Mag_uS_dust' in sagdat.datasetList():
         uMag = sagdat.readDataset('Magnitudes/Mag_uS_dust', idxfilter=massflt)
         rMag = sagdat.readDataset('Magnitudes/Mag_rS_dust', idxfilter=massflt)
      elif 'SED/Magnitudes/Mag_ext_id119_tot_r' in sagdat.datasetList():
         uMag = sagdat.readDataset('SED/Magnitudes/Mag_ext_id119_tot_r')
         rMag = sagdat.readDataset('SED/Magnitudes/Mag_ext_id121_tot_r')
      else:
         print("ERROR: SDSS-u,r total magnitudes not found in output")
         return
      
      urcolor = uMag - rMag
      rMagplot = rMag - 5*np.log10(un.h/0.7)
      del uMag, rMag

      colorlim = divpar[0]-0.2444*np.tanh((rMagplot+divpar[1])/1.09)
      del rMagplot
      
      redflt = urcolor > colorlim
      del colorlim

      logSM  = np.log10(StellarMass*un.h)  # h^-2
      del StellarMass 
  
      Vol = sagdat.boxSizeMpc**3
      bins = np.arange(6, 14, 0.2)

      phi_i, edges = np.histogram(logSM, bins=bins)
      phi = phi_i.astype(float)/Vol/(edges[1]-edges[0])

      phi_i, edges = np.histogram(logSM[redflt], bins=bins)
      phi_red = phi_i.astype(float)/Vol/(edges[1]-edges[0])

      del logSM

   else: #read histograms from file:
      f = h5py.File(readfile, "r")
      phi     = f['RedFrac/phi']
      phi_red = f['RedFrac/phi_red']
      edges   = f['RedFrac/edges']
   
   if savefile:
      f = h5py.File(savefile)
      if 'RedFrac' in f.keys(): del f['RedFrac']
      f['RedFrac/phi']     = phi
      f['RedFrac/phi_red'] = phi_red
      f['RedFrac/edges']   = edges

   # Observations compilation
   fname = datapath+"redfracz0.1_Henriques15Compilation.dat"
   x1, x2, ob_y, ob_e = np.loadtxt(fname, unpack=True)
   ob_x = (x1 + x2)/2.

   # and the plot:
   pl.figure()
   x = (edges[1:]+edges[:-1])/2.0
   l1, = pl.plot(x[phi_red>0], phi_red[phi_red>0]/phi[phi_red>0], '-r', lw=2, label='SAG')
   l2 = pl.errorbar(ob_x, ob_y, yerr=ob_e, fmt='^k', label="Henriques et al (2015), z=0.1")
   pl.xlabel(r'$\log_{10}(M_\star {\rm h}^{-2}[{\rm M}_\odot])$')
   pl.ylabel(r'$\Phi_{\rm red} / \Phi_{\rm total} $')
   pl.xticks([8,9,10,11])
   pl.xlim((7.5,11.5))
   pl.ylim((0,1))
   pl.legend(loc='upper left', frameon=False, handles=[l2,l1])
   pl.tight_layout()
   pl.savefig(outpath+'/RedFraction.eps')
   print("Done!")
   if getPlot: return pl
   pl.clf()
   return 
   
