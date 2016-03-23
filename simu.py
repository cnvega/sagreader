#!/usr/bin/env python

import SAGreader as sag
import SAGplots
import SAGplots_evol
import os

def plots_MDPL():

   outfolder = "plots/SAG-7.86-MDPL"
   outdat = outfolder+"/data.h5"
   
   if not os.path.exists(outfolder):
      os.makedirs(outfolder)

   inpath = "/fast_scratch3/cnvega/MDPL/SAG/SAG7.86-completed/snapshot_125"
   data = sag.SAGcollection(inpath, 1000.0)   


   # 1st: smf
   SAGplots.SMF(data, outfolder, savefile=outdat, redshift=0)

   # 2nd: morph
   SAGplots.FracMorph(data, outfolder, savefile=outdat)

   # BH-B
   SAGplots.BHBulge(data, outfolder, savefile=outdat)

   # SFR density:
   snaplist = sag.SnapList("/fast_scratch2/cnvega/MDPL/snapidzred.txt")
   SAGplots_evol.SFRvol_z(data, snaplist, outfolder, savefile=outdat)

   # High-z SMF:
   #data_tmp = data.select_redshift(0.99, 1.1)
   #SAGplots.SMF(data_tmp, outfolder, savefile=outdat, redshift=1)

   #data_tmp = data.select_redshift(1.99, 2.1)
   #SAGplots.SMF(data_tmp, outfolder, savefile=outdat, redshift=2)

   data.clear()

   sagfolder = "/fast_scratch3/cnvega/MDPL/SAG/SAG7.86-completed/snapshot_092"
   data = sag.SAGcollection(sagfolder, 1000.0) 
   SAGplots.SMF(data, outfolder, savefile=outdat, redshift=1)
   data.clear()

   sagfolder = "/fast_scratch3/cnvega/MDPL/SAG/SAG7.86-completed/snapshot_075"
   data = sag.SAGcollection(sagfolder, 1000.0) 
   SAGplots.SMF(data, outfolder, savefile=outdat, redshift=2)
   data.clear()



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


def plots_MBII():

   outfolder = "plots/SAG-7.92-MBII"
   outdat = outfolder+"/data.h5"
   
   inpath = "/data/MBII/SalidaSAM/SAG-7.92"
   data = sag.SAGcollection(inpath, 100)   

   if not os.path.exists(outfolder):
      os.makedirs(outfolder)

   # 1st: smf
   SAGplots.SMF(data, outfolder, savefile=outdat)

   # 2nd: morph
   SAGplots.FracMorph(data, outfolder, savefile=outdat)

   # BH-B
   SAGplots.BHBulge(data, outfolder, savefile=outdat)

   # SFR density:
   #snaplist = SAGreader.SnapList("/data/MBII/snapidzred.txt")
   #SAGplots_evol.SFRdensity_z(data, snaplist, outfolder, savefile=outdat)

   # High-z SMF:
   #data_tmp = data.select_redshift(0.99, 1.1)
   #SAGplots.SMF(data_tmp, outfolder, savefile=outdat, redshift=1)

   #data_tmp = data.select_redshift(1.99, 2.1)
   #SAGplots.SMF(data_tmp, outfolder, savefile=outdat, redshift=2)

   data.clear()



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
   #plots_stand("SAG-7.c87")
   #plots_MBII()
   plots_MDPL_replot()

