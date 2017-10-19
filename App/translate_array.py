import numpy as np
import logging

logger = logging.getLogger(__name__)

def translate_array(arr,shifts):
    """
    'translate' performs integer translation on a 2D numpy array.

    Params:
      arr - a 2D nparray
      shifts - (xshift, yshift)
    Returns
      new array of same shape and datatype as arr

    Comments:
      Positive values shift up and left
      Areas outside of the original array are filled with np.NAN
        Correction - they would be if we could store NAN values, otherwise they are coerced to 0

    """
    #check array is only 2D
    if not len(arr.shape)==2:
        logger.debug('Only 2D arrays are supported')
        raise RuntimeError

    if not len(shifts)==2:
        logger.debug('Only 2D arrays can be translated')
        raise RuntimeError

    if not shifts[0]==int(shifts[0]) and shifts[1]==int(shifts[1]):
        logger.debug('Only integer translation is supported')
        raise RuntimeError

    new_data = np.zeros(arr.shape, arr.dtype)
    new_data[:]=np.NAN
    s=arr.shape
    new_data[max(0,0-shifts[0]):min(s[0],s[0]-shifts[0]),
             max(0,0-shifts[1]):min(s[1],s[1]-shifts[1])] \
        = arr[max(0,shifts[0]):min(s[0],s[0]+shifts[0]),
              max(0,shifts[1]):min(s[1],s[1]+shifts[1])]

    # Seems we cant store NAN values in an integer array
    if np.any(new_data<0):
        logger.warn('Coercing array values to 0')
        new_data[new_data<0] = 0
    return new_data


if __name__=="__main__":
    orig = np.zeros((5,5),np.uint8)
    orig[2,2] = 1
    new = translate_array(orig, (2,0))
    print(new[0,2] == 1)
