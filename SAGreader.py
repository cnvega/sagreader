#! /usr/bin/env python

import h5py
import numpy as np

class _UnitMass:
   def __init__(self, in_gr):
      self.gr   = in_gr
      self.kg   = in_gr/1e3
      self.Msun = in_gr/1.989e33 
      self.e10Msun = self.Msun/1e10

class _UnitLength:
   def __init__(self, in_cm):
      self.cm = in_cm
      self.m  = in_cm/1e2
      self.kpc = in_cm/3.085678e21
      self.pc  = self.kpc*1e3
      self.Mpc = self.kpc/1e3
      
class _UnitVelocity:
   def __init__(self, in_cm_per_s):
      self.cm_per_s = in_cm_per_s
      self.m_per_s  = in_cm_per_s/1e2
      self.km_per_s = self.m_per_s/1e3
      self.km_per_hr = self.km_per_s*3600.0

class _UnitTime:
   def __init__(self, in_s):
      self.s = in_s
      self.hr = in_s/3600.0
      self.yr = self.hr/24.0/365.0
      self.Myr = self.yr/1e6
      self.Gyr = self.Myr/1e3

class Units:
   """
   The class 'Units' has as attributes all the possible unit conversion in
   a specific SAG data file.
   """
   def __init__(self, length_in_cm, mass_in_gr, vel_in_cm_per_s, h):
      self.mass   = _UnitMass(mass_in_gr)
      self.length = _UnitLength(length_in_cm)
      self.velocity = _UnitVelocity(vel_in_cm_per_s)
      self.time   = _UnitTime(length_in_cm/vel_in_cm_per_s)
      self.h = h

   def unit(self, unit_tag):
      tag = unit_tag.lower()
      if tag in ['msun', 'ms', 'sollarmasses']:
         return self.mass.Msun/self.h
      elif tag in ['msun/h', 'ms/h', 'sollarmasses/h']:
         return self.mass.Msun
      elif tag in ['1e10msun', '1e10ms']:
         return self.mass.e10Msun/self.h
      elif tag in ['1e10msun/h', '1ems/h']:
         return self.mass.e10Msun
      elif tag in ['mpc']:
         return self.length.Mpc/self.h
      elif tag in ['mpc/h']:
         return self.length.kpc
      elif tag in ['kpc']:
         return self.length.kpc/self.h
      elif tag in ['kpc/h']:
         return self.length.kpc
      elif tag in ['km/s']:
         return self.velocity.km_per_s
      elif tag in ['km/hr']:
         return self.velocity.km_per_hr
      elif tag in ['h', 'hubble_h', 'hubble']:
         return self.h
      else:
         raise KeyError


class SAGdata:
   """
   The class 'SAGdata' stores a collection of hdf5 output files
   created by the SAG code. It can extract a particular array from
   all the stored files and returns a unique numpy array with the 
   requested data.
   """
   def __init__(self, simname, boxSizeMpc):
      """
      It creates an empty collection of files.
      """
      self.simname = str(simname)
      self.filenames = []
      self.dataList = []
      self.nfiles = 0
      self.boxSizeMpc = boxSizeMpc
      self.reduced = False


   def clear(self):
      self.simname = ""
      del self.filenames[:]
      self.nfiles = self.boxSizeMpc = 0
      self.reduced = False
      for fsag in self.dataList:
         fsag.close()
        

   def addFile(self, filename):
      """
      It adds an hdf5 file to the object.
      """
      try:
         sag = h5py.File(filename, "r")
      except IOError:
         print("Cannot load file: '"+filename+"'")
         return
      self.filenames.append(filename)
      self.dataList.append(sag)
      self.nfiles += 1
     
      if 1 == self.nfiles:
         try:
            attr = self.dataList[0].attrs['REDUCED_HDF5']
            if attr == 'YES':
               self.reduced = True
         except KeyError:
            pass


   def readDataset(self, dsname, idxfilter=[]):
      """
      It returns a unique numpy array of the requested dataset only
      if it exists in all loaded SAG files.

      The idxfilter can be created with numpy.where(condition), for example:
      >>> types = d.readDataset("Type")
      >>> row, col = numpy.where(types == 0)
      >>> discMass  = d.readDataset("DiscMass", idxfilter=row)
      >>> pos = d.readDataset("Pos", idxfilter=row)
      """
      for i, sag in enumerate(self.dataList):
         dsarr = np.array(sag.get(dsname))
         if None == dsarr.all():
            print("Dataset '"+dsname+"' not present in "+self.filenames[i])
            return None
         if 0 == i:
            nparr = dsarr
         else:
            nparr = np.concatenate([nparr, dsarr])
      if 0 != len(idxfilter):
         tmp = nparr[idxfilter]
         del nparr
         nparr = tmp
            
      return nparr


   def readAttr(self, attname, fnum=0):
      """
      It returns the value of the requested attribute from a particular
      file of the list.
      """
      try:
         attr = self.dataList[fnum].attrs[attname]
         return attr
      except KeyError:
         print("Attribute '"+attname+"' not present in "
               +self.filenames[fnum])
         return None


   def readUnits(self, fnum=0):
      """
      It returns an instance of the 'Units' class, with all the unit conversions
      of the data found in the firts hdf5 file of the list.
      """
      if 0 < self.nfiles:
         if not self.reduced:
            m_in_g = float(self.dataList[fnum].attrs["UnitMass_in_g"])
            l_in_cm = float(self.dataList[fnum].attrs["UnitLength_in_cm"])
            vel_in_cm_s = float(self.dataList[fnum].attrs["UnitVelocity_in_cm_per_s"])
         else:
            m_in_g = 1.989e33  # Msun
            l_in_cm = 3.085678e21  # kpc
            vel_in_cm_s = 1e5      # km/s

         h = float(self.readAttr('Hubble_h'))
         
         units = Units(l_in_cm, m_in_g, vel_in_cm_s, h)
         return units
      else:
         return None

   def datasetList(self, fnum=0, group="/"):
      ks = []
      for tag in self.dataList[fnum][group].keys():
         if type(self.dataList[fnum][group+tag]) is h5py._hl.dataset.Dataset:
            ks.append(group+tag)
         elif type(self.dataList[fnum][group+tag]) is h5py._hl.group.Group:
            tmp = self.datasetList(fnum, group=group+tag+"/")
            ks += tmp
      return ks


   def _gal_idxs(self, ids, dsname):

      if type(ids) != list: ids = [ids]
      idxs = []
      boxes = []
      for i in range(self.nfiles):
         dset = self.dataList[i][dsname]
         tmp = np.where(np.in1d(dset, ids, assume_unique=True))[0]
         idxs += tmp.tolist()
         for _ in range(len(tmp)): boxes.append(i)
      return np.array(idxs), np.array(boxes) 
  

   def getGalaxies(self, dslist='all'):
      if dslist == 'all':
         dslist = self.datasetList()
      gal = {}
      for dstag in dslist:
         if type(self.dataList[0][dstag]) is h5py._hl.dataset.Dataset:
            gal[dstag] = self.readDataset(dstag)
      return gal


   def getGalaxies_by_ids(self, ids, dslist='all'):
      """
      It returns a dictionary with the different datasets for all the requested
      galaxies.
      """
      if dslist == 'all':
         dslist = self.datasetList()
         dslist.remove('Histories/DeltaT_List')
      # retrieve indexes of the galaxies:
      idname = 'GalaxyID' if self.reduced else 'UID'
      idxs, boxes = self._gal_idxs(ids, idname)
      gal = {}
      for dstag in dslist:
         if type(self.dataList[0][dstag]) is h5py._hl.dataset.Dataset:
            dims = self.dataList[0][dstag].shape[1]
            l = np.zeros((len(idxs),dims), dtype=self.dataList[0][dstag].dtype)
            for i in range(self.nfiles):
               l_idx = (boxes == i)
               l[l_idx] = self.dataList[i][dstag][:][idxs[l_idx]]
            gal[dstag] = l
      return gal


class SAGcollection():
   """
   A collection of SAGdata objects of SAG. It assumes the outputs are ordered in
   different directories (one per snapshot/redshift).
   """

   def __init__(self, filename, boxSizeMpc=0):
      """
      It creates a new catalog collection from a specified directory
      """
      self.dataList  = []
      self.snaptag   = []
      self.nfiles    = []
      self.redshift  = []
      self.nz        = 0
      self.boxSizeMpc = 0
      self.zminidx    = -1
     
      if 0 != boxSizeMpc:
         simname = 'SAG_sim'
         self.boxSizeMpc = boxSizeMpc
      else:
         # This file must exist!
         simdat = open(filename+"/simdata.txt")
         simname = simdat.readline()
         self.boxSizeMpc = simdat.readline()
         simdat.close()

      import os
      ls = os.listdir(filename)
      ls.sort()
      for name in ls:
         if name.split("_")[0] == "snapshot":
            snap = name.split("_")[1]
            lss = os.listdir(filename+"/"+name)
            lss.sort()
            filesindir = 0
            for h5name in lss:
               if h5name.split(".")[-1] in ['hdf5','h5'] and \
                  h5name.split("_")[0] == 'gal':
                  # This is a SAG file!:
                  filesindir += 1

                  if filesindir == 1: ## add a new redshift to the collection:
                     self.dataList.append(SAGdata(simname, self.boxSizeMpc))
                     self.snaptag.append(snap)
                     self.nfiles.append(0)
                     self.nz += 1
                  
                  sagname = filename+'/'+name+'/'+h5name
                  print('Opening file '+sagname)
                  self.dataList[self.nz-1].addFile(sagname)
                  self.nfiles[self.nz-1] += 1
      # and the corresponding redshifts: 
      for i in range(self.nz):
         self.redshift.append(float(self.dataList[i].readAttr('Redshift')))


      # If the outputs are not ordered in subfolders, then they are mixed up:
      filesindir = 0
      if 0 == self.nz:
         for name in ls:
            if name.split(".")[-1] in ['hdf5','h5'] and \
               name.split("_")[0] == 'gal':
               
               filesindir += 1
               snap =  name.split("_")[1]
               if snap == 'itf':
                     snap = name.split("_")[2]
               try:
                  idx = self.snaptag.index(snap)
               except ValueError:
                  self.dataList.append(SAGdata(simname, self.boxSizeMpc))
                  self.snaptag.append(snap)
                  self.nfiles.append(0)
                  idx = self.nz
                  self.nz += 1

               sagname = filename+'/'+name
               print('Opening file '+sagname)
               self.dataList[idx].addFile(sagname)
               self.nfiles[idx] += 1
  
               # and the redshift:
               if 1 == self.nfiles[idx]: 
                  self.redshift.append(float(self.dataList[idx].readAttr('Redshift')))

      self.zminidx = self.redshift.index(min(self.redshift))
      self.reduced = self.dataList[0].reduced


   def clear(self):
      del self.snaptag[:]
      del self.nfiles[:]
      del self.redshift[:]
      self.nz = 0
      self.boxSizeMpc = 0
      self.zminidx = -1
      for sagd in self.dataList:
         sagd.clear()
   

   def _lookup_z(self, zlow, zhigh):
      """
      It returns a list with the redshifts of the collection that are in the
      range zlow <= Z <= zhigh.
      """
      zl = []
      for z in self.redshift:
         if zlow <= z <= zhigh:
            zl.append(z)
      return zl


   def readDataset(self, dsname, multiSnaps=False, zrange=None, **kwargs):
      """
      It searches for an unique or a set of redshifts or boxes and returns the 
      requested datasets.
      """
      for key in kwargs.keys():
         if key not in ['idxfilter']:
            raise KeyError(key)
      
      if multiSnaps: print("Warning: requesting multiple snaps!")

      if not zrange:
         if not multiSnaps:
            # the lowest redshift (hopefully z=0):
            iarr = [self.zminidx]
         else:
            iarr = [i for i in range(self.nz)]
      else:
         zl = zrange[0]
         zh = zrange[1]
         if multiSnaps:
            iarr = []
            for z in self._lookup_z(zl, zh):
               iarr.append(self.redshift.index(z))
         else:
            # search for the lowest match:
            zarr = self._lookup_z(zl, zh)
            iarr = [self.redshift.index(min(zarr))]
     
      if 'idxfilter' in kwargs.keys():
         flt = kwargs['idxfilter']
      else:
         flt = []
      # Now we have the list of redshifts we are going to use, let's concatenate the
      # datasets:
      for k, i in enumerate(iarr):
         dsarr = self.dataList[i].readDataset(dsname)
         if 0 == k:
            nparr = dsarr
         else:
            nparr = np.concatenate([nparr, dsarr])

      if 0 != len(flt):
         tmp = nparr[flt]
         del nparr
         nparr = tmp

      return nparr
        

   # All the files should have the same attributes and units, so these return
   # the ones found in the first file at the lowest redshift.
   def readAttr(self, attname):
      """
      It returns a requested attribute.
      """
      return self.dataList[self.zminidx].readAttr(attname)
   
   def readUnits(self):
      """ 
      It return an 'Units' object with the unit conversions of the catalog.
      """
      return self.dataList[self.zminidx].readUnits() 

   def datasetList(self):
      """
      It returns the dataset list of the files, following the groups recursively.
      """
      return self.dataList[self.zminidx].datasetList()


   def getGalaxies_by_ids(self, ids, dslist='all', multiSnaps=False, zrange=None):
      """
      It returns a dictionary with the different datasets for all the requested
      galaxies.
      """
      if multiSnaps: print("Warning: requesting multiple snaps!")

      if not zrange:
         if not multiSnaps:
            # the lowest redshift (hopefully z=0):
            iarr = [self.zminidx]
         else:
            iarr = [i for i in range(self.nz)]
      else:
         zl = zrange[0]
         zh = zrange[1]
         if multiSnaps:
            iarr = []
            for z in self._lookup_z(zl, zh):
               iarr.append(self.redshift.index(z))
         else:
            # search for the lowest match:
            zarr = self._lookup_z(zl, zh)
            iarr = [self.redshift.index(min(zarr))]

      if dslist == 'all':
         dslist = self.dataList[iarr[0]].datasetList()
         if 'Histories/DeltaT_List' in dslist:
            dslist.remove('Histories/DeltaT_List')

      # here we need to verify the datasets are present in all the snaps.
      if multiSnaps:
         if len(iarr) > 1:
            ds1 = set(dslist)
            ds2 = set(self.dataList[iarr[1]].datasetList())
            dslist = list(ds1.intersection(ds2))
            
      # now we collect the galaxies:

      gal = { 'ngal': np.array([0]) }
      for i in iarr:
         gtmp = self.dataList[i].getGalaxies_by_ids(ids, dslist)
         if 0 == gal['ngal']:
            gal.update(gtmp)
            gal['ngal'] += len(gal[dslist[0]])
            if len(iarr) > 1:
               gal['redshift'] = np.repeat(self.redshift[i],len(gal[dslist[0]]))
         else:
            for ds in dslist:
               gal[ds] = np.concatenate([gal[ds], gtmp[ds]])
            gal['ngal'] += len(gtmp[dslist[0]])
            tmp = np.repeat(self.redshift[i],len(gtmp[dslist[0]]))
            gal['redshift'] = np.concatenate([gal['redshift'], tmp])

      return gal
      
      

class SAGreader():
   """
   Input SAG Mock
   """
   def __init__(self, filename=''):
      self.catalog = SAGcollection(filename)
      
