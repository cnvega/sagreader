#!/usr/bin/env python

import SAGreader as sag
import SAGplots
import SAGplots_evol
import os

def plots_MDPL():

   outfolder = "plots/SAG-7.128"
   outdat = outfolder+"/data.h5"
   
   if not os.path.exists(outfolder):
      os.makedirs(outfolder)

   inpath = "/data2/MDPL/SalidaSAM/SAG-7.128/snapshots/snapshot_125"
   data = sag.SAGcollection(inpath, boxSizeMpc=1000., keepOpen=False)   
   
   SAGplots.set_style("book", Wfrac=0.75)
   
   ## 1st: smf
   SAGplots.SMF(data, outfolder, savefile=outdat, redshift=0)

   ## 2nd: morph
   #SAGplots.FracMorph(data, outfolder, savefile=outdat)

   ## BH-B
   #SAGplots.BHBulge(data, outfolder, savefile=outdat)
   #
   ## Tully-Fisher
   ##SAGplots.TullyFisher(data, outfolder, savefile=outdat)

   # High-z SMF:
   #data_tmp = data.select_redshift(0.99, 1.1)
   #SAGplots.SMF(data_tmp, outfolder, savefile=outdat, redshift=1)

   #data_tmp = data.select_redshift(1.99, 2.1)
   #SAGplots.SMF(data_tmp, outfolder, savefile=outdat, redshift=2)


   # SFR density:

   #sagfolder = "/fast_scratch3/cnvega/MDPL/SAG/SAG7.86-completed"
   #data = sag.SAGcollection(sagfolder, 1000.0) 

   #snaplist = sag.SnapList("/fast_scratch2/cnvega/MDPL/snapidzred.txt")

   #SAGplots_evol.SFRvol_z(data, snaplist, outfolder, savefile=outdat)

   #data_tmp = data.select_redshift(0.14, 0.155)  
   #SAGplots.GasFrac(data_tmp, outfolder, savefile=outdat)
   #SAGplots.SFRF(data, outfolder, savefile=outdat)
   #SAGplots.MstarMhalo(data, outfolder, savefile=outdat)
   
   data.clear()



def plots_stand(modelid):
   
   outfolder = "plots/Stand/"+modelid
   outdat = outfolder+"/data.h5"

   if not os.path.exists(outfolder):
      os.makedirs(outfolder)

   salidaSAM = "/data/Stand/SalidaSAM/"+modelid+"/Subgalaxies"

   data = sag.SAGcollection(salidaSAM, boxSizeMpc=150.0, keepOpen=False)

   # 1st: smfs
   SAGplots.SMF(data, outfolder, savefile=outdat)

   #data_tmp = data.select_redshift(0.98, 1.05)
   #SAGplots.SMF(data_tmp, outfolder, savefile=outdat, redshift=1)

   #data_tmp = data.select_redshift(1.98, 2.1)
   #SAGplots.SMF(data_tmp, outfolder, savefile=outdat, redshift=2)

   #data_tmp = data.select_redshift(2.98, 3.1)
   #SAGplots.SMF(data_tmp, outfolder, savefile=outdat, redshift=3)

   # 2nd: morph
   SAGplots.FracMorph(data, outfolder, savefile=outdat)

   # BH-B
   SAGplots.BHBulge(data, outfolder, savefile=outdat)

   # T-F
   SAGplots.TullyFisher(data, outfolder, savefile=outdat)
   
   # SFR density:

   snaplist = sag.SnapList("/data/Stand/outputs_STAND.dat")

   SAGplots_evol.SFRvol_z(data, snaplist, outfolder, savefile=outdat)

   data.clear()

def replot_folder(outfolder):

   outdat = outfolder+"/data.h5"
   
   if not os.path.exists(outfolder):
      print(outfolder+" does not exist!")
      return 

   #SAGplots.set_style('book', Wfrac=0.75)
   SAGplots.set_style("mnras", Hratio=0.75)

   ### Full set of plots:
   #SAGplots.SMF(None, outfolder, readfile=outdat, redshift=0)
   #SAGplots.SMF(None, outfolder, readfile=outdat, redshift=1)
   #SAGplots.SMF(None, outfolder, readfile=outdat, redshift=2)
   #SAGplots.FracMorph(None, outfolder, readfile=outdat)
   #SAGplots.BHBulge(None, outfolder, readfile=outdat)
   #SAGplots.TullyFisher(None, outfolder, readfile=outdat)
   #SAGplots.CMD(None, outfolder, readfile=outdat, divpar=(1.8, 18.7))
   #SAGplots.RedFraction(None, outfolder, readfile=outdat, divpar=(1.8, 18.7))
   #SAGplots.GasFrac(None, outfolder, readfile=outdat)
  
   #SAGplots.SFRF(None, outfolder, readfile=outdat)
   
   SAGplots.MstarMhalo(None, outfolder, readfile=outdat)
   SAGplots.MstarMhalo(None, outfolder, readfile=outdat, galTypes=[0])
   SAGplots.MstarMhalo(None, outfolder, readfile=outdat, galTypes=[1])


   #SAGplots_evol.SFRvol_z(None, None, outfolder, readfile=outdat)

if __name__ == '__main__':
   #replot_folder("plots/SAG-7.128")
   plots_MDPL()
   #plots_stand("SAG-7.96-c05abr16")
   #plots_MDPL_replot()

