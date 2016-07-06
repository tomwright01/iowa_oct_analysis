import numpy as np
import logging
from libtiff import TIFF

logger = logging.getLogger(__name__)

def readTiffFileStack(fname):
    """Read a tiff image stack (i.e. exported from imageJ)
    and return a ndarray"""
    im = TIFF.open(fname)
    width = im.getfield('ImageWidth')
    height = im.getfield('ImageLength')
    img_description = im.getfield('ImageDescription')
    img_description = dict(x.split('=') for x in img_description.split('\n') if x)
    frames = int(img_description['slices'])
    
    data_array = np.zeros((frames, height, width))
    for idx,img in enumerate(im.iter_images()):
        data_array[idx,:,:] = img
        
    info = {'ascans' : width,
            'bscans' : frames}
    
        
    return {'headers' : info,
            'image_data' : data_array}

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    fname = '../Data/Sample/Bioptigen/Oct Intensity Stack.tif'
    result = readTiffFileStack(fname)    
    