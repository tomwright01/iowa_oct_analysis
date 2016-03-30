import numpy as np
import re
import os
import datetime
import logging

logger = logging.getLogger(__name__)

def readCirrusOct(fname):
    """Read a Cirrus .img file
        Params:
        fname - full path to .img file

        Returns:
        {'headers':{},
        'image_data':np.array(nBscans,scan_depth,nAscans),}

        Examples
        headers, image_data = readCirrusOct(fname)"""
    
    if not fname.lower().endswith('.img'):
        raise RuntimeError('Only cirrus .img files are supported')
    
    # parse the filename to get the scan size (voxels)
    info = parseCirrusFilename(fname)
    
    image = np.fromfile(fname, dtype='uint8')

    expected_size = info['ascans'] * info['bscans'] * 1024
    assert len(image) == expected_size, 'Raw file size doesnt match layers size'
    
    image = image.reshape((info['bscans'],
                           1024,
                           info['ascans']))

    # flip array vertically, needed fo compatibility with layer data
    image = image[:,::-1,:]
    # and flip horizontally
    image = image[:,:,::-1]        
    
    return {'headers':info,
            'image_data':image}

def parseCirrusFilename(fname):
    """Parse a cirrus .img filename
    
    Params:
    fname - filename (can include path)
    
    Returns:
    """
    dirname = os.path.dirname(fname)
    filename = os.path.basename(fname)
    
    match_str = '(\d+)x(\d+)_(\d{2}-\d{2}-\d{4})_(\d{2}-\d{2}-\d{2})_(OS|OD)_(sn\d+)'
    parts = re.search(match_str, filename)
    
    if parts is None:
        logger.warning('Failed to parse filename: %s', fname)
        return None
    dateParts = parts.group(3).split('-')
    timeParts = parts.group(4).split('-')
    
    testDateTime = datetime.datetime(int(dateParts[2]),
                                     int(dateParts[0]),
                                     int(dateParts[1]),
                                     int(timeParts[0]),
                                     int(timeParts[1]),
                                     int(timeParts[2]))
    
    results = {'dirname':dirname,
              'filename':filename,
              'ascans':int(parts.group(1)),
              'bscans':int(parts.group(2)),
              'testTimeDate':testDateTime,
              'eye':parts.group(5),
              'serial':parts.group(6)}
    
    return results
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    fname = '../Data/Sample/Cirrus/Sample1/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z.img'
    result = readCirrusOct(fname)