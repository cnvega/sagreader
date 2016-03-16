#! /usr/bin/env python

import sys
import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pylab as pl
import numpy as np
import math
import h5py

from matplotlib import rcParams

#rcParams['font.size'] = 16

def SMF(sagdat, outpath, savefile=None, readfile=False, redshift=0):

   print("### Stellar Mass Function")

   Hubble_h = sagdat.readAttr('Hubble_h')
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

      else :
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
      StellarMass /= Hubble_h
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
      if 'SMF' in f.keys(): del f['SMF']
      f['SMF/z'+str(z)+'/bins'] = bins
      f['SMF/z'+str(z)+'/phi']   = phi
      f.close()

   # Finally, the plot:
   pl.figure()
   x = (bins[1:]+bins[:-1])/2.0
   pl.plot(x[phi>0], np.log10(phi[phi>0]), '-k', lw=2, label='SAG')
   pl.errorbar(ob_x, np.log10(ob_y) , yerr=[ob_ed, ob_eu], fmt='^b',
              label="Henriques et al (2013), $z=0$")
   pl.xlabel(r'$\log_{10}(M_{\rm stellar}[{\rm M}_\odot])$', fontsize=16)
   pl.ylabel(r'$\log_{10}(\Phi[{\rm h}^{-3} {\rm Mpc}^{-3} /\log_{10} M])$', 
             fontsize=16)
   pl.xlim((7,12.5))
   pl.legend(loc=3, frameon=False, fontsize=12)
   #pl.axis('scaled')
   pl.tight_layout()
   pl.savefig(outpath+'/SMF_z'+str(z)+'.eps')
   pl.clf()
   print("Done!") 

def FracMorph(sagdat, outpath, savefile=None, readfile=False,
               threshSE=0.85, threshIS=0.0, SEDmagfilter='NONE'):
   
   print("### Morphological Fraction")

   nbins = 40
   mrange = [7,13]
   
   datapath = './Data/'
   Hubble_h = sagdat.readAttr('Hubble_h')

   if not readfile:
      
      # First we check if there are an available Magnitude for pre-selecting the 
      # galaxies:

      if 'Magnitudes/Mag_V_dust1' in sagdat.datasetList():
         MagV = sagdat.readDataset('Magnitudes/Mag_V_dust1')
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
      f.close()

   if savefile:
      f = h5py.File(savefile)
      if 'Morph' in f.keys(): del f['Morph']
      f['Morph/total_counts'] = tot
      f['Morph/bins'] = bin_tot
      f['Morph/irr'] = irr 
      f['Morph/sp']  = sp  
      f['Morph/el']  = el  
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
   pl.plot(x, f_sp, "b--", label="Spirals", linewidth=3)
   pl.plot(x, f_el, "r-.", label="Ellipticals", linewidth=3)
   pl.plot(x, f_irr, "k:", label="Irregulars", linewidth=3)
   pl.legend(frameon=False, fontsize=12, loc='upper left')

   pl.errorbar(obs_sp[0]+np.log10(Hubble_h), obs_sp[1], yerr=obs_sp[2], fmt="db")
   pl.errorbar(obs_el[0]+np.log10(Hubble_h), obs_el[1], yerr=obs_el[2], fmt="or")
   pl.errorbar(obs_irr[0]+np.log10(Hubble_h), obs_irr[1], yerr=obs_irr[2], fmt="^k")

   pl.xlabel(r"$\log_{10}(M_{\rm stellar} [h^{-1} {\rm M}_\odot])$", fontsize=16)
   pl.ylabel("Fraction", fontsize=16)
   pl.xlim((8.,12.))
   pl.ylim((-0.05, 1.05))
   pl.tight_layout()
   pl.savefig(outpath+'/FracMorph.eps')
   pl.clf()
   print("Done!")


def BHBulge(sagdat, outpath, savefile=None, readfile=False, SEDmagfilter='NONE'):

   print("### Black hole vs Bulge mass relationship")

   from matplotlib import colors, ticker, cm
   from numpy import ma

   Hubble_h = sagdat.readAttr('Hubble_h')
   datapath = './Data/'

   xr = [8,13]
   yr = [5,11]

   if not readfile:
      if 'Magnitudes/Mag_V_dust1' in sagdat.datasetList():
         MagV = sagdat.readDataset('Magnitudes/Mag_V_dust1')
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
      BulgeMass *= un.mass.Msun/Hubble_h
      BHMass *= un.mass.Msun/Hubble_h
      
      nz = (BulgeMass>0)&(BHMass>0)

      H, xedges, yedges = np.histogram2d(np.log10(BulgeMass[nz]), np.log10(BHMass[nz]), 
                                         range=[xr,yr], bins=80)
      del BulgeMass, BHMass

   else: #read histograms from file:
      f = h5py.File(readfile, "r")
      H = f['BHBulge/histogram']
      xedges = f['BHBulge/Xedges_bulge']
      yedges = f['BHBulge/Yedges_bh']
   
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
               fmt="ob", label="Kormendi & Ho (2013)")
   ax.errorbar(obM[0], obM[2], xerr=obM[1], yerr=[obM[4], obM[3]], 
               fmt="^g", label="McConnell & Ma (2013)")

   # This is just for not including the error bars in the legend
   handles, labels = ax.get_legend_handles_labels()
   handles = [h[0] for h in handles]
   ax.legend(handles, labels, frameon=False, fontsize=12, loc="upper left",
             numpoints=1)

   #cs = ax.contourf(x, y, H, vmin=1, locator=ticker.LogLocator(), cmap=cm.Reds)
   cs = ax.imshow(H.T, origin="low", cmap=cm.YlOrBr, norm=colors.LogNorm(vmin=1),
                  extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]], 
                  aspect='auto', interpolation='nearest')
   fig.colorbar(cs)


   ax.set_xlabel(r"$\log_{10}(M_{\rm Bulge} [{\rm M}_\odot])$", fontsize=16)
   ax.set_ylabel(r"$\log_{10}(M_{\rm BH} [{\rm M}_\odot])$", fontsize=16)
   
   ax.set_xlim(xr)
   ax.set_ylim(yr)

   pl.tight_layout()

   fig.savefig(outpath+"/BHBulge.eps")
   pl.clf()
   print("Done!")
   

def TullyFisher(sagdat, outpath, savefile=None, readfile=False, SEDmagfilter='NONE'):

   print("### Tully Fisher relationship")
   from matplotlib import colors, ticker, cm
   from numpy import ma

   Hubble_h = sagdat.readAttr('Hubble_h')
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
         Mag_rS = sagdat.readDataset('Magnitudes/Mag_rS_dust1', idxfilter=flt)
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
   ax.legend(handles, labels, frameon=False, fontsize=12, loc="upper left",
             numpoints=1)

   #cs = ax.contourf(x, y, H, vmin=1, locator=ticker.LogLocator(), cmap=cm.Greens)
   cs = ax.imshow(H.T, origin='lower', 
                  cmap=cm.Greens, norm=colors.LogNorm(vmin=1),
                  extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]], 
                  aspect='auto', interpolation='nearest')
   fig.colorbar(cs)

   #ax.plot(np.log10(Vrot), Magr, ".k", markersize=0.1)

   ax.set_xlabel(r"$\log_{10}(V_{\rm flat} [{\rm km/s}])$", fontsize=16)
   ax.set_ylabel(r"$M_{\rm r} - 5 \log_{10} h$", fontsize=16)
   
   ax.set_xlim(xr)
   ax.set_ylim([yr[1],yr[0]])

   pl.tight_layout()

   fig.savefig(outpath+"/TullyFisher.eps")
   pl.clf()
   print("Done!")
   
   
