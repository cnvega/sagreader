#! /usr/bin/env python
# coding: utf-8

## @file SAGreader.py
## @author Cristian A. Vega Martínez <cnvega(at)fcaglp.unlp.edu.ar>
##
## @brief SAG output reader module.
##
## A collection of classes to load and extract data from standard and
## reduced h5df output files of the SAG code for galaxy formation and 
## evolution. 


import h5py
import numpy as np

class _UnitMass:
   """
   Internal class which stores the most used mass conversion constants.
   """
   def __init__(self, in_gr):
      self.gr   = in_gr
      self.kg   = in_gr/1e3
      self.Msun = in_gr/1.989e33 
      self.e10Msun = self.Msun/1e10

class _UnitLength:
   """
   Internal class which stores the most used length conversion constants.
   """
   def __init__(self, in_cm):
      self.cm = in_cm
      self.m  = in_cm/1e2
      self.kpc = in_cm/3.085678e21
      self.pc  = self.kpc*1e3
      self.Mpc = self.kpc/1e3
      
class _UnitVelocity:
   """
   Internal class which stores the most used velocity conversion constants.
   """
   def __init__(self, in_cm_per_s):
      self.cm_per_s = in_cm_per_s
      self.m_per_s  = in_cm_per_s/1e2
      self.km_per_s = self.m_per_s/1e3
      self.km_per_hr = self.km_per_s*3600.0

class _UnitTime:
   """
   Internal class which stores the most used time conversion constants.
   """
   def __init__(self, in_s):
      self.s = in_s
      self.hr = in_s/3600.0
      self.yr = self.hr/24.0/365.0
      self.Myr = self.yr/1e6
      self.Gyr = self.Myr/1e3

class Units:
   """ Unit conversion constants.

The class 'Units' has as attributes all the possible unit conversion constants
in a specific SAG data file or collection. The object can be used as conversion
factor of the desired unit without considering the 'h' scaling, or applying a full
conversion by using the unit() method.

Example:

un = Unit(l_cm, m_gr, vel_cm_s)

StellarMass_in_gr = stellarMass*un.mass.gr

Positions_in_mpc   = positions*un.length.Mpc
   """

   def __init__(self, length_in_cm, mass_in_gr, vel_in_cm_per_s, h):
      """ Constructor.
      
This four parameters are needed.

@param length_in_cm: length conversion constant in cm
@param mass_in_gr: mass conversion constant in gr
@param vel_in_cm_per_s: velocity conversion constant in cm_s
@param h: Hubble constant in units of 100 km_s_Mpc
      """

      self.mass   = _UnitMass(mass_in_gr)
      self.length = _UnitLength(length_in_cm)
      self.velocity = _UnitVelocity(vel_in_cm_per_s)
      self.time   = _UnitTime(length_in_cm/vel_in_cm_per_s)
      self.h = h

   def unit(self, unit_tag):
      """ Get conversion factor from tag.

It calculates and returns the conversion constant of a requested unit 
by considering the 'h' scaling.

@param unit_tag string with the desired unit (examples: 'Msun', 'kpc/h').
A KeyError is raised when the tag is not found.

@return The requested conversion factor.
      """
      tag = unit_tag.lower()
      if tag in ['msun', 'ms', 'sollarmasses']:
         return self.mass.Msun/self.h
      elif tag in ['msun/h', 'ms/h', 'sollarmasses/h']:
         return self.mass.Msun
      elif tag in ['1e10msun', '1e10ms']:
         return self.mass.e10Msun/self.h
      elif tag in ['1e10msun/h', '1e10ms/h']:
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
   """ One snapshot SAG output containing multiple files/boxes.

The class 'SAGdata' stores a collection of hdf5 output files
created by the SAG code. It can extract a particular array from
all the stored files and returns a unique numpy array with the 
requested data. All the hdf5 files have to be added manually.
"""

   def __init__(self, simname, boxSizeMpc, keepOpen=True):
      """ Constructor.

It creates an empty collection of files. The parameter values are not used in
this class, but are useful when plotting the data.
      
@param simname String containing a tag name for the simulation (unused).
@param boxSizeMpc Box size of the simulation box (unused).
@param keepOpen (optional) It controls if the hdf5 files are always opened or only
when needed.

      """
      ## Reference name of the simulation.
      self.simname = str(simname)
      ## List with the names of opened hdf5 files.
      self.filenames = []
      ## List with the opened hdf5 files of the SAG output.
      self.dataList = []
      ## Number of loaded files.
      self.nfiles = 0
      ## Box size of the simulation in Mpc.
      self.boxSizeMpc = boxSizeMpc
      ## Tag indicating if the hdf5 files are in the reduced format or not.
      self.reduced = False
      ## Tag indicating if the hdf5 file are opened only when needed or not.
      self.keepOpen = keepOpen


   def clear(self):
      """ Object cleaner.

It deletes all the internal lists and it closes all the loaded files.
      """
      self.simname = ""
      del self.filenames[:]
      self.nfiles = self.boxSizeMpc = 0
      self.reduced = False
      self.keepopen = True
      if not self.keepOpen:
         for fsag in self.dataList:
            fsag.close()
        

   def addFile(self, filename):
      """ Add one file to the object.

It adds and opens one hdf5 file to the internal lists. An Input/Output error is
raised if the file cannot be loaded.

@param filename Name of the file to be loaded. The path can be absolute or
relative.
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
            if type(attr) is np.bytes_:
               attr = attr.decode()
            if 'YES' == attr:
               self.reduced = True
         except KeyError:
            pass
      
      if not self.keepOpen:
         for i in range(self.nfiles):
            self.dataList[i].close()


   def readDataset(self, dsname, idxfilter=[]):
      """ Dataset reader.

It returns a unique numpy array of the requested dataset only
if it exists in all loaded SAG files. It creates a concatenated array
considering all the files.

@param dsname String with name of the desired dataset.
@param idxfilter (optional) numpy array for filtering the requested data. 
It can be created with numpy.where(condition), for example:

> types = d.readDataset("Type")

> row, col = numpy.where(types == 0)

> discMass  = d.readDataset("DiscMass", idxfilter=row)

> pos = d.readDataset("Pos", idxfilter=row)

@return A numpy array with the requested dataset.
      """
      for i, fname in enumerate(self.filenames):
         if self.keepOpen: sag = self.dataList[i]
         else: sag = h5py.File(self.filenames[i], "r")
         dsarr = np.array(sag.get(dsname))
         if None == dsarr.all():
            print("Dataset '"+dsname+"' not present in "+fname)
            return None
         if 0 == i:
            nparr = dsarr
         else:
            nparr = np.concatenate([nparr, dsarr])
         if not self.keepOpen:
            sag.close()
      if 0 != len(idxfilter):
         tmp = nparr[idxfilter]
         del nparr
         nparr = tmp
            
      return nparr


   def readAttr(self, attname, fnum=0):
      """ Attribute reader.

It returns the value of the requested attribute from a particular
file of the list. A KeyError is raised if not found.

@param attname String with the name of the desired attribute. 
@param fnum (optional) File number, starting from zero, from which
the attribute is extracted. The first file is used as default.

@return Content of the requested attribute.
      """
      if self.keepOpen: sag = self.dataList[fnum]
      else: sag = h5py.File(self.filenames[fnum], "r")
      try:
         attr = sag.attrs[attname]
      except KeyError:
         print("Attribute '"+attname+"' not present in "
               +self.filenames[fnum])
         attr = None
      if not self.keepOpen: sag.close()
      return attr


   def readUnits(self, fnum=0):
      """ Unit reader.

It returns an instance of the 'Units' class, with all the unit conversions
of the data.

@param fnum (optional) File number, starting from zero, from which
the units are extracted. The first file is used as default.

@return A 'Unit' object with all the conversion constant of the loaded data.
      """
      if self.keepOpen: sag = self.dataList[fnum]
      else: sag = h5py.File(self.filenames[fnum], "r")
      if 0 < self.nfiles:
         if not self.reduced:
            m_in_g = float(sag.attrs["UnitMass_in_g"])
            l_in_cm = float(sag.attrs["UnitLength_in_cm"])
            vel_in_cm_s = float(sag[fnum].attrs["UnitVelocity_in_cm_per_s"])
         else:
            m_in_g = 1.989e33  # Msun
            l_in_cm = 3.085678e21  # kpc
            vel_in_cm_s = 1e5      # km/s

         h = float(self.readAttr('Hubble_h'))
         
         units = Units(l_in_cm, m_in_g, vel_in_cm_s, h)
      else: 
         units = None
      
      if not self.keepOpen: sag.close()
      return units


   def datasetList(self, fnum=0, group="/"):
      """ Get the list of datasets.

It recursively extracts the list of datasets from a group of a particular loaded
file.

@param fnum (optional) File number, starting from zero, from which
the dataset keys are extracted. The first file is used as default.
@param group (optional) Group in the hdf5 hierarchy. The group is used as default.

@return A List with all the names of the datasets found. 
      """
      if self.keepOpen: sag = self.dataList[fnum]
      else: sag = h5py.File(self.filenames[fnum], "r")
      ks = []
      for tag in sag[group].keys():
         if type(sag[group+tag]) is h5py._hl.dataset.Dataset:
            ks.append(group+tag)
         elif type(sag[group+tag]) is h5py._hl.group.Group:
            tmp = self.datasetList(fnum, group=group+tag+"/")
            ks += tmp
      if not self.keepOpen: sag.close()
      return ks


   def _gal_idxs(self, ids, dsname):
      """ Get index of galaxies from their IDs.

It returns the internal location of the galaxies identified by a list of ids

@param ids A unique or a list of galaxy IDs to search.
@param dsname The name of the dataset in which the galaxy ids are stored.

@return (idxs, boxes) It returns two numpy arrays specifying the location of the
galaxies. An 'indexes' array with the location of the galaxies in their corresponding
files, and a 'boxes' array indicating the file to which each galaxy belongs. 
      """
      if self.keepOpen: sag = self.dataList[fnum]
      else: sag = h5py.File(self.filenames[fnum], "r")
      if type(ids) != list: ids = [ids]
      idxs = []
      boxes = []
      for i in range(self.nfiles):
         dset = sag[dsname]
         tmp = np.where(np.in1d(dset, ids, assume_unique=True))[0]
         idxs += tmp.tolist()
         for _ in range(len(tmp)): boxes.append(i)
      if not self.keepOpen: sag.close()
      return np.array(idxs), np.array(boxes) 
  

   def getGalaxies(self, dslist='all'):
      """ Retrieve galaxies.

It creates a dictionary with all the requested datasets including all the galaxies
found in the loaded files. 

@param dslist (optional) A string or a list of strings with the names of the datasets
to be included. The default behavior is to include all the datasets found in the
files.

@return A dictionary with one numpy array for each requested dataset. The names of
the datasets are preserved as tags of the dictionary.
      """
      if dslist == 'all':
         dslist = self.datasetList()
      gal = {}
      for dstag in dslist:
         #if type(self.dataList[0][dstag]) is h5py._hl.dataset.Dataset:
         gal[dstag] = self.readDataset(dstag)
      return gal


   def getGalaxies_by_ids(self, ids, dslist='all'):
      """ Retrieve galaxies by their ids.

It creates a dictionary with all the requested datasets including ONLY the galaxies
with the desired ids.

WARNING: This method could be time-consuming.

@param ids A unique or a list of galaxy IDs to search.
@param dslist (optional) A string or a list of strings with the names of the datasets
to be included. The default behavior is to include all the datasets found in the
files.

@return A dictionary with one numpy array for each requested dataset. The names of
the datasets are preserved as tags of the dictionary.
      """
      if dslist == 'all':
         dslist = self.datasetList()
         if not self.reduced:
            dslist.remove('Histories/DeltaT_List')
      # retrieve indexes of the galaxies:
      idname = 'GalaxyID' if self.reduced else 'UID'
      idxs, boxes = self._gal_idxs(ids, idname)
      gal = {}
      if self.keepOpen: sag = self.dataList[0]
      else: sag = h5py.File(self.filename[0], "r")
      for dstag in dslist:
         if type(sag[dstag]) is h5py._hl.dataset.Dataset:
            dims = sag[dstag].shape[1]
            l = np.zeros((len(idxs),dims), dtype=sag[dstag].dtype)
            for i in range(self.nfiles):
               l_idx = (boxes == i)
               l[l_idx] = sag[dstag][:][idxs[l_idx]]
            gal[dstag] = l
      if not self.keepOpen: sag.close()
      return gal


class SAGcollection():
   """ A collection of SAG outputs containing multiple boxes and redshifts. 

The class 'SAGcollection' is a compound of 'SAGdata' objects classified by
their redshifts. It behaves as a 'SAGdata' object but included the redshift
management of the data, including the dataset and galaxy requests. The 'multisnap'
flag included as argument in most of the method controls if the outputs must extract
and return the data from different redshifts or just the lowest one. A z range can be
specified in each case.
   """

   def __init__(self, filename, boxSizeMpc=0, keepOpen=True):
      """ It creates a new catalog collection.

It reads a specific directory and creates a collection of SAG outputs classified and
ordered by their redshifts. 

@param filename Input directory in which the SAG files are looked and loaded.
Those hdf5 files can share the same directory in which case the 
classification and ordering is made by
reading the 'Redshift' attribute of each one, or can be separated by
redshift/snapshot in
different folders named 'snapshot_NNN', where NNN is the snapshot number.
@param boxSizeMpc (optional) The box size of the simulation in Mpc. If this argument
in omitted, it is loaded by default a file called 'simdata.txt' 
from the catalogs folder. The file must contain just two lines: a tag name of the
simulation and its box size in Mpc.
@param keepOpen (optional) It controls if the hdf5 files are always opened or only
when needed.
      """

      ## List with the 'SAGdata' objects for each redshift/snapshot.
      self.dataList  = []
      ## List with the tag name of each snapshot.
      self.snaptag   = []
      ## List with the number of files loaded in each snapshot.
      self.nfiles    = []
      ## List with the corresponding redshift of each snapshot.
      self.redshift  = []
      ## Number of redshits/snapshots loaded.
      self.nz        = 0
      ## Box size of the simulation in Mpc.
      self.boxSizeMpc = 0
      ## Index in which the lowest redsfhit is located.
      self.zminidx    = -1
      ## Tag indicating if the hdf5 file are opened only when needed or not.
      self.keepOpen = keepOpen
     
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
                     self.dataList.append(SAGdata(simname, self.boxSizeMpc,
                                                  self.keepOpen))
                     self.snaptag.append(snap)
                     self.nfiles.append(0)
                     self.nz += 1
                  
                  sagname = filename+'/'+name+'/'+h5name
                  print('Loading file '+sagname)
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
                  self.dataList.append(SAGdata(simname, self.boxSizeMpc,
                                               self.keepOpen))
                  self.snaptag.append(snap)
                  self.nfiles.append(0)
                  idx = self.nz
                  self.nz += 1

               sagname = filename+'/'+name
               print('Loading file '+sagname)
               self.dataList[idx].addFile(sagname)
               self.nfiles[idx] += 1
  
               # and the redshift:
               if 1 == self.nfiles[idx]: 
                  self.redshift.append(float(self.dataList[idx].readAttr('Redshift')))

      self.zminidx = self.redshift.index(min(self.redshift))
      self.reduced = self.dataList[0].reduced


   def clear(self):
      """ Object cleaner.

It deletes all the internal lists and it closes all the loaded files.
      """
      del self.snaptag[:]
      del self.nfiles[:]
      del self.redshift[:]
      self.nz = 0
      self.boxSizeMpc = 0
      self.zminidx = -1
      for sagd in self.dataList:
         sagd.clear()
   

   def _lookup_z(self, zlow, zhigh):
      """ Search a redshift range.

It searches for the loaded redshifts which are in the specified range.

@param zlow Lower redshift of the range.
@param zhigh Higher redshift of the range.
@return A list with the redshifts of the collection that are in the
range zlow <= Z <= zhigh.
      """
      zl = []
      for z in self.redshift:
         if zlow <= z <= zhigh:
            zl.append(z)
      return zl


   def select_redshift(zmin, zmax):
      zm = min(self._lookup_z(zmin, zmax))
      idx = self.redshift.index(zm)
      return self.dataList[idx]



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
      
      
class SnapList():
   def __init__(self, fname):
      f = open(fname, "r")
      self.snap = []
      scales = []
      redshifts = []
      snp = 0
      for line in f.readlines():
         if line[0] == '#': continue
         # check the format of the list:
         values = line.split() 
         if 5 == len(values):
            snap, a, z = values[0:3]
         if 1 == len(values):
            snap = snp
            snp += 1
            a = values[0]
            z = 1./float(a) - 1.
         self.snap.append(int(snap))
         scales.append(float(a))
         redshifts.append(float(z))
      self.scales = np.array(scales)
      self.redshifts = np.array(redshifts)
   def __getitem__(self, k):
      idx = self.snap.index(k)
      return self.scales[idx], self.redshifts[idx]
       

class SAGreader():
   """
   Input SAG Mock
   """
   def __init__(self, filename=''):
      self.catalog = SAGcollection(filename)
      
