import numpy as np
import logging

logger = logging.getLogger(__name__)

def readBioptigenOct(fname):
    """Read a biotigen .oct file
    Params:
      fname - full path to .oct file
      
    Returns:
      np.array
    """
    # create a dict to store the file header info
    # each item has the datatype and the length, 
    # if length is 0 then the value of data_length is used

    header_defs = {'FRAMEHEADER':[np.uint32,1],
                   'FRAMECOUNT':[np.uint32,1],
                   'LINECOUNT':[np.uint32,1],
                   'LINELENGTH':[np.uint32,1],
                   'SAMPLEFORMAT':[np.uint32,1],
                   'DESCRIPTION':[np.uint8,0],
                   'XMIN':[np.float64,1],
                   'XMAX':[np.float64,1],
                   'XCAPTION':[np.uint8,0],
                   'YMIN':[np.float64,1],
                   'YMAX':[np.float64,1],
                   'YCAPTION':[np.uint8,0],
                   'SCANTYPE':[np.uint32,1],
                   'SCANDEPTH':[np.float64,1],
                   'SCANLENGTH':[np.float64,1],
                   'AZSCANLENGTH':[np.float64,1],
                   'ELSCANLENGTH':[np.float64,1],
                   'OBJECTDISTANCE':[np.float64,1],
                   'SCANANGLE':[np.float64,1],
                   'SCANS':[np.uint32,1],
                   'FRAMES':[np.uint32,1],
                   'FRAMESPERVOLUME':[np.uint32,1],
                   'DOPPLERFLAG':[np.uint32,1],
                   'CONFIG':[np.uint8,0]}
    headers = {}
    
    with open(fname, 'rb') as fid:
        magic_number = np.fromfile(fid, np.uint16, 2)
        version_number = np.fromfile(fid, np.uint16, 1)
        key_length = np.fromfile(fid, np.uint32, 1)
        key = ''.join([chr(v) for v in np.fromfile(fid, np.uint8, key_length)])
        data_length = np.fromfile(fid, np.uint32, 1)
        
        if key != 'FRAMEHEADER':
            logger.error('Failed loading frame header for file %s'.format(fname))
            raise RuntimeError('Error loading frame header')
        
        header_flag = False # set to True when all header keys are read
        while not header_flag:
            key_length = np.fromfile(fid, np.uint32, 1)
            key = ''.join([chr(v) for v in np.fromfile(fid, np.uint8, key_length)])
            data_length = np.fromfile(fid, np.uint32, 1)
            
            if key in header_defs.keys():
                if header_defs[key][1]:
                    headers[key] = np.fromfile(fid, header_defs[key][0], 1)
                else:
                    headers[key] = np.fromfile(fid, header_defs[key][0], data_length)
            else:
                header_flag = True

    pass
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    fname = '../Data/Sample/Bioptigen/sample_OD_V_12x12_0_0000036.OCT'
    readBioptigenOct(fname)
    