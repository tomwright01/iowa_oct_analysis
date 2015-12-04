from nose.tools import *
import numpy as np
import App.OctLayers as Oct

def test_center_data_1():
    oct = Oct.OctLayers()
    oct.data = np.ones((1,5,5))
    oct.center_x = 1
    oct.center_y = 1
    oct.centerData()
    
    target = np.empty((1,5,5))
    target[:] = np.NAN
    target[0,1:5,1:5] = np.ones((4,4))
    assert np.all(oct.data==target)
    
def test_center_data_2():
    oct = Oct.OctLayers()
    oct.data = np.ones((1,5,5))
    oct.center_x = 3
    oct.center_y = 3
    oct.centerData()
    
    target = np.empty((1,5,5))
    target[:] = np.NAN
    target[0,0:4,0:4] = np.ones((4,4))
    assert np.all(oct.data==target)
    
if __name__=='__main__':
    logging.basicConfig(level=logging.DEBUG)
    test_center_data()