#!/usr/bin/env python

import SAGreader as sag
import SAGplots
import SAGplots_evol
import os
   
import h5py
import numpy as np
import matplotlib.mlab as mlab
import math


def compare_models():
   
   outfolder = "plots/paper/"

   if not os.path.exists(outfolder):
      print(outfolder+" does not exist!")
      return 

   outdat = "plots/SAG-7.128-paper/data.h5"
   SAGplots.set_style("mnras", Hratio=0.75)

   SAGplots.ExtraDatSMF = True

   # this file was created with the whole z=0 snapshot of the \beta_F=1.3 model: 
   f = h5py.File("SAG_aF1.3.h5", "r")

   ## -------------------------------------------------
   ## SMF at z=0
   p = SAGplots.SMF(None, outfolder, readfile=outdat, redshift=0, getPlot=True)
   ax = p.gca()  
   ax.lines[-2].set_label(r'SAG, $\beta_\mathrm{F}=1.99$')
   
   bins = f['SMF/z0/bins'][:]
   phi  = f['SMF/z0/phi'][:]
   
   x = (bins[1:]+bins[:-1])/2.0
   p.plot(x[phi>0], np.log10(phi[phi>0]), '--b', lw=1.5, 
          label=r'SAG$_{\beta1.3}$', zorder=9)
   p.legend(loc=3, frameon=False, fontsize='small', handlelength=2,
            labelspacing=0.3)
   p.savefig(outfolder+'/SMF_z0_comp.eps')

   ## -------------------------------------------------
   ### BH-Bulge mass relationship:
   p = SAGplots.BHBulge(None, outfolder, readfile=outdat, getPlot=True)

   H = f['BHBulge/histogram'][:]
   xedges = f['BHBulge/Xedges_bulge'][:]
   yedges = f['BHBulge/Yedges_bh'][:]

   H = np.ma.masked_where(H<=0, H)
   x = (xedges[1:]+xedges[:-1])/2
   y = (yedges[1:]+yedges[:-1])/2

   X,Y = np.meshgrid(x,y)
   
   newH = H.T/(xedges[1]-xedges[0])/(yedges[1]-yedges[0])
   ticks = np.array([0.01, 0.19, 0.26, 0.38, 0.68, 0.95, 0.997])*newH.max()

   p.text(0.8, 0.1, r'$z = 0$', transform=p.gca().transAxes)

   p.contour(X, Y, newH, levels=ticks, zorder=11, colors='b', linewidths=0.5,
             linestyles=':')
   p.tight_layout()
   p.savefig('BHBulge_comp.eps')

   ## -------------------------------------------------
   ### Gas fraction of the galaxies:
   p = SAGplots.GasFrac(None, outfolder, readfile=outdat, getPlot=True)
  
   H = f['GasFrac/histogram'][:]
   xedges = f['GasFrac/Xedges_mstar'][:]
   yedges = f['GasFrac/Yedges_gasfrac'][:]
   sag_x    = f['GasFrac/sag_x'][:]
   sag_mean = f['GasFrac/sag_mean'][:]
   sag_std  = f['GasFrac/sag_std'][:]
   
   H = np.ma.masked_where(H<=0, H)
   x = (xedges[1:]+xedges[:-1])/2
   y = (yedges[1:]+yedges[:-1])/2
   
   X,Y = np.meshgrid(x,y)
   newH = H.T/(xedges[1]-xedges[0])/(yedges[1]-yedges[0])
   ticks = np.array([0.01, 0.19, 0.26, 0.38, 0.68, 0.95, 0.997])*newH.max()

   p.contour(X, Y, newH, levels=ticks, zorder=9, colors='b', linewidths=0.5,
             linestyles=':')
   
   p.text(0.8, 0.8, r'$z = 0$', transform=p.gca().transAxes)

   p.plot(sag_x, sag_mean, 'b--', lw=1.5, label='SAG')
   #p.plot(sag_x, sag_mean+sag_std, 'b--', lw=0.5)
   #p.plot(sag_x, sag_mean-sag_std, 'b--', lw=0.5)
   
   p.tight_layout()
   p.savefig('GasFrac_comp.eps')
   
   f.close()


if __name__ == '__main__':
   compare_models()

