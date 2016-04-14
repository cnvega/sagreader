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
   """ SAG output containing multiple files/boxes of one snapshot.

The class 'SAGdata' stores a collection of hdf5 output files
created by the SAG code. It can extract a particular array from
all the stored files and it returns an unique numpy array with the 
requested data. All the hdf5 files have to be added manually.
"""

   def __init__(self, simname, boxSizeMpc, keepOpen=True):
      """ Constructor.

It creates an empty collection of files. The required parameter values
are not used in this class, but are useful when plotting the data.
      
@param simname String containing a tag name for the simulation (unused).

@param boxSizeMpc Box size of the simulation box (unused).

@param keepOpen (optional) It controls the accessibility of the hdf5
files by controlling if they are always opened or only when needed.

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
The constructor arguments are preserved.
      """
      #self.simname = ""
      del self.filenames[:]
      self.nfiles = 0
      #self.boxSizeMpc = 0
      self.reduced = False
      if self.keepOpen:
         for fsag in self.dataList:
            fsag.close()
      #self.keepopen = True
        

   def addFile(self, filename):
      """ Add one file to the object.

It adds one hdf5 file to the internal lists. An Input/Output error is
raised if the file cannot be opened.

@param filename Name of the file to be loaded. The path can be absolute or
relative.
      """
      try:
         sag = h5py.File(filename, "r")
      except IOError:
         print("Cannot load file: '"+filename+"'")
         return
      self.filenames.append(filename)
      self.nfiles += 1

      if 1 == self.nfiles:
         try:
            attr = sag.attrs['REDUCED_HDF5']
            if type(attr) is np.bytes_:
               attr = attr.decode()
            if 'YES' == attr:
               self.reduced = True
         except KeyError:
            pass
      
      if self.keepOpen:
         self.dataList.append(sag)
      else:
         sag.close()

   def readDataset(self, dsname, idxfilter=[]):
      """ Dataset reader.

It returns a unique numpy array of the requested dataset only
if it exists in all loaded SAG files. It creates a concatenated array
considering all the files.

@param dsname String with name of the desired dataset.

@param idxfilter (optional) A numpy array for filtering the requested 
data in the first index of the array. It can be created with
numpy.where(condition), for example:

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

It returns the content of the requested attribute from a particular
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
         if isinstance(sag[group+tag], h5py._hl.dataset.Dataset):
            ks.append(group+tag)
         elif isinstance(sag[group+tag], h5py._hl.group.Group):
            tmp = self.datasetList(fnum, group=group+tag+"/")
            ks += tmp
      if not self.keepOpen: sag.close()
      #deleting the initial '/':
      for i in range(len(ks)):
         if ks[i][0] == '/': ks[i] = ks[i][1:]
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
      if type(ids) != list: ids = [ids]
      idxs = []
      boxes = []
      for i in range(self.nfiles):
         if self.keepOpen: sag = self.dataList[i]
         else: sag = h5py.File(self.filenames[i], "r")
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
         if dstag[0] == '/': dstag = dstag[1:]
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
      idname = 'GalaxyStaticID' if self.reduced else 'UID'
      idxs, boxes = self._gal_idxs(ids, idname)
      gal = {}
      if self.keepOpen: s0 = self.dataList[0]
      else: s0 = h5py.File(self.filenames[0], "r")
      for dstag in dslist:
         dims = s0[dstag].shape[1]
         l = np.zeros((len(idxs),dims), dtype=s0[dstag].dtype)
         for i in range(self.nfiles):
            if self.keepOpen: sag = self.dataList[i]
            else: sag = h5py.File(self.filenames[i], "r")
            l_idx = (boxes == i)
            if len(idxs[l_idx] > 0):
               l[l_idx] = sag[dstag][:][idxs[l_idx]]
            if not self.keepOpen: sag.close()
         if dstag[0] == '/': dstag = dstag[1:]
         gal[dstag] = l
      if not self.keepOpen: s0.close()
      return gal


class SAGcollection():
   """ A collection of SAG outputs containing multiple boxes and redshifts. 

The class 'SAGcollection' is a compound of 'SAGdata' objects classified by
their redshifts. It behaves as a 'SAGdata' object but it considers the redshift
management of the data, including the dataset and galaxy requests. The 'multisnap'
flag included as argument in most of the methods controls if the outputs must extract
and return the data from different redshifts or just the lowest one. A redshift 
range can be specified in each case.
   """

   def __init__(self, filename, boxSizeMpc=0, keepOpen=True):
      """ It creates a new catalog collection.

It reads a specific directory and creates a collection of SAG outputs classified and
ordered by their redshifts. 

@param filename Input directory in which the SAG files are looked and loaded.
Those hdf5 files can share the same directory in which case the 
classification and ordering is made by reading the 'Redshift' 
attribute of each one, or can be separated by redshift/snapshot in
different folders named 'snapshot_NNN', where NNN is the snapshot number.

@param boxSizeMpc (optional) The box size of the simulation in Mpc. If this argument
in omitted, a file called 'simdata.txt' is loaded by default 
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
      """ Search for a redshift range.

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


   def select_redshift(self, zmin, zmax):
      """ Get SAGdata with specific redshift.

It creates a reference to the 'SAGdata' object stored in the collection which
redshift matches within the specified range. If more than one redshift is found 
in the range, the lowest is chosen.

@param zmin Minimum redshift of the search range.

@param zmax Maximum redshift of the search range.

@return A 'SAGdata' reference with the requested data.
      """
      zm = min(self._lookup_z(zmin, zmax))
      idx = self.redshift.index(zm)
      return self.dataList[idx]



   def readDataset(self, dsname, multiSnaps=False, zrange=None, **kwargs):
      """ Dataset reader

It returns a unique numpy array of the requested dataset only
if it exists in all the corresponding SAG files, creating a concatenated 
array.

@param dsname String with name of the desired dataset.

@param multiSnaps (optional) It controls if multiple redshift snapshots are
considered to be included in the array.

@param zrange (optional) Redshift range from which the data is extracted. If
multiSnaps=False and more than one redshift is found in the range, 
the lowest redshift is chosen.

@param kwargs Set of optinal parameters of the SAGdata.readDataset method.

@return A numpy array with the requested dataset.
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
      """ Attribute reader.

It returns the content of the requested attribute, which is extracted from the 
first file of the SAGdata object with lowest redshift. 
A KeyError is raised if not found.

@param attname String with the name of the desired attribute. 

@return Content of the requested attribute.
      """
      return self.dataList[self.zminidx].readAttr(attname)
  

   def readUnits(self):
      """ Unit reader.

It returns an instance of the 'Units' class, with all the unit conversions
of the data. It is extracted from the first file of the SAGdata object with 
lowest redshift. 

@return A 'Unit' object with all the conversion constant of the catalog.
      """
      return self.dataList[self.zminidx].readUnits() 


   def datasetList(self):
      """ Get the list of datasets.

It recursively extracts the list of datasets found in the first 
file of the SAGdata object with lowest redshift.

@return A List with all the names of the datasets found. 
      """
      return self.dataList[self.zminidx].datasetList()


   def getGalaxies_by_ids(self, ids, dslist='all', multiSnaps=False, zrange=None):
      """ Retrieve galaxies by their ids.

It creates a dictionary with all the requested datasets including ONLY the galaxies
with the desired ids. A 'redshift' array is included in the dictionary if multiple
snapshots are requested.

WARNING: This method could be time-consuming.

@param ids A unique or a list of galaxy IDs to search.

@param dslist (optional) A string or a list of strings with the names of the datasets
to be included. The default behavior is to include all the datasets found in the
files.

@param multiSnaps (optional) It controls if multiple redshift snapshots are
considered.

@param zrange (optional) Redshift range from which the data is extracted. If
multiSnaps=False and more than one redshift is found in the range, 
the lowest redshift is chosen.

@return A dictionary with one numpy array for each requested dataset. The names of
the datasets are preserved as tags of the dictionary.
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
            gal['ngal'] += len(gtmp[dslist[0]])
            if len(iarr) > 1:
               gal['redshift'] = np.repeat(self.redshift[i],len(gal[dslist[0]]))
         else:
            for ds in dslist:
               gal[ds] = np.concatenate([gal[ds], gtmp[ds]])
            gal['ngal'] += len(gtmp[dslist[0]])
            tmp = np.repeat(self.redshift[i],len(gtmp[dslist[0]]))
            gal['redshift'] = np.concatenate([gal['redshift'], tmp])
         
      gal['redshift'] = gal['redshift'].reshape((len(gal['redshift']),1))
      return gal
      
      
class SnapList():
   """ Scales list

A class for loading, ordering and managing a scales (and redshift) list,
according to their corresponding snapshot numbers.
   """
   def __init__(self, fname):
      """ Constructor

It reads a scales/redshifts list from a file, which can be in any of the 
two most used formats: i) a single column with the scale values in ascending 
order, one line per snapshot; or ii) a five column file in which the first three
columns are the snapshot number, the scale and the redshift (the last two are
omitted here).

@param fname Name of the file for loading the list.

@return A SnapList object with the loaded list.
      """
      ## List of snapshot numbers.
      self.snap = []

      scales = []
      redshifts = []

      f = open(fname, "r")
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

      ## List of scale factors.
      self.scales = np.array(scales)
      ## List of redshifts.
      self.redshifts = np.array(redshifts)


   def __getitem__(self, k):
      """ Get item 

The squared bracket operator receives the snapshot number as argument and return a
tuple with the corresponding scale and redshift.

@param k Snapshot number

@return (a, z) 
      """
      idx = self.snap.index(k)
      return self.scales[idx], self.redshifts[idx]
       

class SAGreader():
   """ Input SAG Mock 

Dummy class used as reference for other external catalogs.
   """
   def __init__(self, filename=''):
      self.catalog = SAGcollection(filename)
      
