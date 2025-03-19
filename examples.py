import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import App.OctLayers as oct

def layerMaskExample(obj):
    """Demonstrates how to get a mask for an OCT layer and apply it to the OCT images
    Note surfaces are numbered from 0 not 1"""
    layer_mask = obj.getOctLayerMask(3,4,True)
    myoct = np.copy(obj.octdata)
    myoct[layer_mask] = 255
    plt.imshow(myoct[0,:,:])
    plt.show()
    
def etdrsThicknessExample(obj):
    """Demonstrates how to get layer thickness values for the 9 etdrs regions
    Note surfaces are numbered from 0 not 1"""
    thickness1 = obj.getEtdrsThickness(0,3,True)
    thickness2 = obj.getEtdrsThickness(3,5,True)
    thickness = pd.concat([thickness1, thickness2], axis=1)
    
    # TODO: make this a proper test
    # layers 0-3 should be thicker than 3-5 (except at the fovea)
    print (thickness.loc['Fovea'][0] < thickness.loc['Fovea'][1])
    print (np.all((thickness.loc[:,0] > thickness.loc[:,1])[1:]))
    print (thickness)


def etdrsIntensityExample(obj):
    """Demonstrates how to get layer intensity values for the 9 etdrs regions
    Note surfaces are numbered from 0 not 1"""
    intensities1 = obj.getEtdrsIntensity(5,6,True)
    intensities2 = obj.getEtdrsIntensity(7,10,True)
    intensities = pd.concat([intensities1, intensities2], axis=1)
    
    # TODO: make this a proper test
    # layers 5 -6 should be dimmer than 7 - 10
    print (np.all(intensities.loc[:,0] < intensities.loc[:,1]))
    print (intensities)
    


if __name__=='__main__':
    plt.gray()
    
    fname = 'Data/Sample/Cirrus/Sample1/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z.img'
    fname_layers = 'Data/Sample/Cirrus/Sample1/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z_Surfaces_Iowa.xml'
    fname_centers = 'Data/Sample/Cirrus/Sample1/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z_GridCenter_Iowa.xml'
    
    data_oct = oct.OctLayers(filename=fname_layers,
                             center_filename=fname_centers,
                             raw_filename=fname)
    data_oct.findFovea()
    data_oct.centerData()
    
    # layerMaskExample(data_oct)
    etdrsThicknessExample(data_oct)
    etdrsIntensityExample(data_oct)