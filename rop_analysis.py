import logging
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import re
import datetime
import fnmatch
import dill
import glob
import collections
import App.OctLayers as Oct
from PIL import Image

if __name__=='__main__':
    logging.basicConfig(level=logging.INFO)

    data_dir = '/mnt/AO/OCT ROP'

    # define the layers of interest
    #    [topSurface,bottomSurface,Name,applyMask]
    #    surfaces are numbered from 0 (not 1 as in segmentation software)
    #layer_set = [(0,10,'Total',False),
                 #(0,1,'RNFL',True),
                 #(1,2,'GCL',True),
                 #(2,3,'IPL',True),
                 #(1,3,'GCL + IPL',True),
                 #(1,4,'GCL + IPL + INL',True),
                 #(0,3,'RNFL + GCL + IPL',True),
                 #(0,4,'RNFL + GCL + IPL + INL',True),
                 #(5,8,'Outer',True),
                 #(6,8,'PhR',True),
                 #(3,4,'INL',True),
                 #(4,5,'OPL',True),
                 #(5,6,'ONL',True),
                 #(1,5,'Inner Complex',True),
                 #(5,8,'Outer Complex',True),
                 #(4,8,'Outer Complex + OPL',True),
                 #(1,8,'Total GCL-Bruchs',True)]
    layer_set = [(0,10,'Total',False),
                 (5,8,'Outer',True),
                 (1,5,'Inner Complex',True),
                 (5,8,'Outer Complex',True),
                 (1,8,'Total GCL-Bruchs',True)]
    # find the analysed files
    file_search = '*_Surfaces_Iowa.xml'
    date_search = re.compile('.*ExamDate_(\d{2}.\d{2}.\d{4}).*')
    # create a named object to hold the file details
    dat = collections.namedtuple('FileInfo','patid, DoT, Eye, Size, Sample, Path, Filename')
    files = []
    for root, dirnames, filenames in os.walk(data_dir):
        for filename in fnmatch.filter(filenames,file_search):
            parts = filename.split('_')
            if parts[3] == '12x12':
                #only append 12 by 12 files
                DateOfTest = re.match(date_search,root)
                if DateOfTest:
                    DateOfTest = datetime.datetime.strptime(DateOfTest.group(1),
                                                            '%m.%d.%Y')
                    DateOfTest = DateOfTest.strftime('%Y-%m-%d')
                files.append(dat._make((parts[0], # patid
                                        DateOfTest, # date of test
                                        parts[1], # eye
                                        parts[3], # size
                                        parts[5], # sample
                                        root, # containing folder
                                        filename)))

    #  setup a data structure to hold the output
    output = {}    
    for layer in layer_set:
        output[layer[2]]=pd.DataFrame()
        
    iCount = 0
    for file in files:
        iCount = iCount + 1
        logging.info('Processing file {} of {}'.format(iCount,len(files)))
        #load the data
        center_filename = '_'.join([file.patid, 
                                    file.Eye, 
                                    'V', 
                                    file.Size, 
                                    '0', 
                                    file.Sample, 
                                    'GridCenter_Iowa.xml'])
        if not os.path.isfile(os.path.join(file.Path, center_filename)):
            logging.warning('Center file not found for patient {} eye {}'.format(file[0],file[1]))
            data = Oct.OctLayers(filename=os.path.join(file.Path,file.Filename),
                                 system='bioptigen',
                                 laterality=file.Eye)
        else:
            data = Oct.OctLayers(filename=os.path.join(file.Path,file.Filename),
                                 system='bioptigen',
                                 laterality=file.Eye,
                                 center_filename= os.path.join(file.Path,center_filename))
            data.centerData()
        
        # convert all to RE orientation
        data.setOrient('OD')
        # process the layer data
        for layer in layer_set:
            etdrs_values = data.getEtdrsThickness(layer[0],
                                                  layer[1],
                                                  True) #ETDRS values as pd.series
            etdrs_values['ID'] = file.patid
            etdrs_values['Eye'] = file.Eye
            etdrs_values['DoT'] = file.DoT
            etdrs_values['Sample'] = int(file.Sample)
            output[layer[2]] = output[layer[2]].append(etdrs_values,ignore_index=True)

    # save the output
    for key in output.keys():
        output[key].to_csv(os.path.join(data_dir,'Output_' + key + '.csv'))
        