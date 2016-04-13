from nose.tools import *
import numpy as np
import App.translate_array as translate
import logging

def test_translate_array():
    orig = np.zeros((5,5),dtype=np.uint8)
    orig[2,2] = 1
    assert(translate.translate_array(orig,(2,0))[0,2]==1)
    assert(translate.translate_array(orig,(-2,0))[4,2]==1)
    assert(translate.translate_array(orig,(0,2))[2,0]==1)
    assert(translate.translate_array(orig,(0,-2))[2,4]==1)
    
    
if __name__=='__main__':
    logging.basicConfig(level=logging.DEBUG)
    test_translate_array()    