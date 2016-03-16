#!/usr/bin/env python

import SAGreader as sag
import SAGplots
import os

def plots_simu():

   outfolder = "plots/MDPL-samp"
   outdat = outfolder+"/data.h5"
   
   data = sag.SAGdata("MDPL", 1000)   
   infile = "gal_125_SAG7.86_BOX_001.hdf5"
   data.addFile(infile)

   if not os.path.exists(outfolder):
      os.makedirs(outfolder)

   # 1st: smf
   SAGplots.SMF(data, outfolder, savefile=outdat)

   # 2nd: morph
   SAGplots.FracMorph(data, outfolder, savefile=outdat)

   # BH-B
   SAGplots.BHBulge(data, outfolder, savefile=outdat)

def plots_stand(modelid):
   
   nboxes = 64
   zsnap = {0:"099"}
   
   outfolder = "plots/"+modelid
   outdat = outfolder+"/data.h5"

   if not os.path.exists(outfolder):
      os.makedirs(outfolder)

   salidaSAM = "/data/Stand/SalidaSAM/"+modelid+"/Subgalaxies"

   #data = sag.SAGdata("Stand", 150.0)

   #for box in range(1,nboxes+1):
   #   infile = salidaSAM+"gal_itf_"+zsnap[0]+"_"+modelid+"_BOX_"+str(box).zfill(3)+".hdf5"
   #   data.addFile(infile)

   data = sag.SAGcollection(salidaSAM, boxSizeMpc=150.0)

   # 1st: smf
   SAGplots.SMF(data, outfolder, savefile=outdat)

   # 2nd: morph
   SAGplots.FracMorph(data, outfolder, savefile=outdat)

   # BH-B
   SAGplots.BHBulge(data, outfolder, savefile=outdat)

   # T-F
   SAGplots.TullyFisher(data, outfolder, savefile=outdat)

   data.clear()


if __name__ == '__main__':
   #plots_simu()
   plots_stand("SAG-7.c87")

