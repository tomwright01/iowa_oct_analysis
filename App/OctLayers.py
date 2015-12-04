import numpy as np
import logging
import os
import fnmatch
import re
import xml.etree.ElementTree

logger = logging.getLogger(__name__)

class OctLayers(object):
    """Class to hold the surface information in xml files from the Iowa reference algorithms
    Usage : oct = OctLayers(filename = <path to surface file>,
                            [center_filename = <path to Grid_Center file>],
                            [laterality = 'OD' | 'OS'])
    If laterality is not provided we will try to autodetect it from the xml file or from the filename
    """
    surface_labels = ['ILM','RNFL-GCL','GCL-IPL','IPL-INL','INL-OPL',
                      'OPL-HFL','BMEIS','IS/OSJ','IB_OPR','IB_RPE','OB_RPE']
    
    def __init__(self, data=None, *args, **kargs):
        if not data:
            self.data = np.empty((11,128,512))
        else:
            self.data = np.array(data)
        
        if 'filename' in kargs.keys():
            self.filename = kargs['filename']
            self.loadXml()
            
        if 'center_filename' in kargs.keys():
            self.center_filename = kargs['center_filename']
            self.loadCenter()
        else:
            self.center_x = None
            self.center_y = None
            
        # can override the autodetected laterality here
        if 'laterality' in kargs.keys():
            self.laterality = kargs['laterality']

    # delgation function
    def __getattr__(self,name):
        """Delegate to NumPy array"""
        try:
            return getattr(self.data, name)
        except:
            raise AttributeError(
                "'Array' object has no attribute {}".format(name))
    
    # delgation function
    def __getitem__(self,key):
        return self.data.__getitem__(key)
    
    # delgation function
    def __setitem__(self,key,value):
        self.data.__setitem__(key,value)    
            
    def loadXml(self):
        """Load the xml file defined in self.filename"""
        if not self.filename:
            logger.error('Source filename not set')
            raise RuntimeError('Filename not set')
        logger.debug('Loading file:{}'.format(self.filename))
        xml_root = xml.etree.ElementTree.parse(self.filename).getroot()
        p = re.compile('.*\((.*)\)')
        for surface in xml_root.findall('surface'):
            # identify which surface this is and assign an index
            # can't use the label in the xml file as these are not contiguous
            surface_name = surface.find('name').text
            logger.debug('Loading surface:{}'.format(surface_name))
            surface_idx = np.NaN
            match = re.match(p, surface_name)
            if match:
                try:
                    surface_idx = self.surface_labels.index(match.group(1))
                except ValueError:
                    logger.debug('Surface {} not defined'.format(match.group(1)))
                    break
            else:
                logger.warning('Failed to identify surface:{}'.format(surface_name))
                break
            logger.debug('Surface index:{}'.format(surface_idx))
            # loop through all the bscans
            surface_bscans = surface.findall('bscan')
            for bscan_idx in range(128):
                bscan = surface_bscans[bscan_idx]
                self.data[surface_idx,bscan_idx,:] = \
                    [int(z.text) for z in bscan.findall('z')]
        
        undef_mask = np.zeros(self.data.shape,dtype=np.bool)
        undef_xml = xml_root.find('undefined_region')
        
        if undef_xml is not None:
            for ascan in undef_xml.findall('ascan'):
                x = int(ascan.find('x').text)
                y = int(ascan.find('y').text)
                undef_mask[:,y,x] = True
        self.data = np.ma.MaskedArray(self.data, mask = undef_mask)

        laterality = xml_root.find('scan_characteristics').find('laterality').text
        if laterality.upper in ['OD','OS']:
            self.laterality = laterality.upper
        else:
            # not defined in the xml file, see if we can extract from the filename
            p = re.compile('(OD|OS)')
            m = re.search(p,self.filename)
            if m:
                self.laterality = m.group(0)
                
    def loadCenter(self):
        """Load the xml file defined in self.center_filename"""
        if not self.center_filename:
            logger.error('Center filename not set')
            raise RuntimeError('Filename not set')
        
        xml_root = xml.etree.ElementTree.parse(self.center_filename)
        c = xml_root.find('center')
        self.center_x = int(c.find('x').text)
        self.center_y = int(c.find('y').text)
        
    def centerData(self):
        """data matrix is shifted so that defined coordinates are in the center of the array
        Out of bounds areas are marked as np.NAN and are masked"""
        
        assert self.center_x and self.center_y, 'Center coords not set'
        centered_data = np.ma.MaskedArray(np.empty(self.data.shape))
        centered_data[:] = np.NAN
        row_shift = self.center_y - int(self.data.shape[1] / 2)
        col_shift = self.center_x - int(self.data.shape[2] / 2)
        
        if row_shift <= 0:
            # shift data down
            t_idx = 0
            b_idx = self.data.shape[1]  - abs(row_shift)
            target_t_idx = abs(row_shift)
        else:
            # shift data up
            t_idx = row_shift
            b_idx = self.data.shape[1]
            target_t_idx = 0
            
        if col_shift <= 0:
            # shift data right
            l_idx = 0
            r_idx = self.data.shape[2] - abs(col_shift)
            target_l_idx = abs(col_shift)
        else:
            # shift data right
            l_idx = col_shift
            r_idx = self.data.shape[2]
            target_l_idx = 0
            
        valid_data = self.data[:,t_idx:b_idx, l_idx:r_idx]
        centered_data[:,
                      target_t_idx:(target_t_idx + valid_data.shape[1]),
                      target_l_idx:(target_l_idx + valid_data.shape[2])] = \
            valid_data
        centered_data.mask = centered_data.mask + np.isnan(centered_data) 
        self.data = centered_data
    
    def setOrient(self,target='OD'):
        assert self.laterality, 'Original laterality not defined'
        if self.laterality != target:
            self.data = np.fliplr(self.data)
            self.laterality = target
    
    def getThickness(self,surface1, surface2):
        return np.abs(self.data[surface1,:,:] - self.data[surface2,:,:])
    
class OctCollection(object):
    
    def __init__(self,folder = None, laterality='OD'):
        assert laterality in ['OD','OS'], 'Laterality mus be OS or OD'
        self.data = []
        self.laterality = laterality
        self.srcFolder = folder
        if folder:
            self.loadFolder()
        
    def loadFolder(self):
        assert self.srcFolder, 'Source folder not set'
        surface_files = []
        center_files = []
        for dirname, _, filenames in os.walk(self.srcFolder):
            for name in fnmatch.filter(filenames,'*Surfaces_Iowa.xml'):
                surface_files.append((dirname,name))
            for name in fnmatch.filter(filenames,'*GridCenter_Iowa.xml'):
                center_files.append((dirname,name))                
        
        p = re.compile('(.*)_Surfaces_Iowa.xml')
        for surface_file in surface_files:
            basename = re.match(p, surface_file[1]).group(1)
            try:
                i = [i for i,v in enumerate(center_files) if re.search(basename,v[1])][0]
            except IndexError:
                i = None
            
            if i is not None:
                data = OctLayers(filename = os.path.join(surface_file[0], surface_file[1]),
                                 center_filename = os.path.join(center_files[i][0],center_files[i][1]))
                logger.debug('Centering data')
                data.centerData()
            else:
                logger.warning('GridCenter file not found for surface file:{}'.format(surface_file[1]))
                data = OctLayers(filename = os.path.join(surface_file[0], surface_file[1]))
                
            data.setOrient(self.laterality)
            self.data.append(data)
            
    def getLayerValues(self, surface1, surface2, value_type='mean'):
        assert self.data, 'No data loaded'
        assert value_type in ['mean','stdev'], 'Invalid value_type requested'
        layers = np.ma.MaskedArray(np.empty((len(self.data),
                                             self.data[0].shape[1],
                                             self.data[1].shape[2])))
        for idx in range(len(self.data)):
            layer = self.data[idx].getThickness(surface1, surface2)
            layers[idx,:,:] = layer
        
        if value_type == 'mean':
            return layers.mean(axis=0)
        
        if value_type == 'stdev':
            return layers.std(axis=0)